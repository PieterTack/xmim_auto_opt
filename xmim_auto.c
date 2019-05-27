#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <xraylib.h>
#include <glib.h>
#include <glib/gstdio.h>
#include <xmi_msim.h>
#include <float.h>

#define XMI_PYMCA_START_CONC 0.0001
#define XMI_PYMCA_MAX_ITERATIONS 100
#define XMO_CONV_THRESHOLD 0.05

struct xmo_spe_data {
	int nchannels;
	double *mca_data;
};

int xmo_read_spe(char* spe_file, struct xmo_spe_data *spe) {
	FILE *filePtr;
	char buffer[1024];
	unsigned long int nchannels[2];
	int i;

	if ( (filePtr = fopen(spe_file,"r")) == NULL){
		fprintf(stderr,"Could not open file %s\n", spe_file);
		return -1;
	}

	do {
		if (fgets(buffer,1024,filePtr) == NULL) {
			fprintf(stderr,"ERROR: An error occured while reading %s\n",spe_file);
			return -1;
		}
	} while (strstr(buffer,"$DATA:") == NULL);
	fscanf(filePtr,"%lu %lu",&(nchannels[0]),&(nchannels[1]));

	spe->nchannels = (int) (nchannels[1]-nchannels[0]+1);
	spe->mca_data = (double *) malloc(sizeof(double)*spe->nchannels);
	for (i = 0; i < spe->nchannels; i++) {
		if(fscanf(filePtr,"%lf",&(spe->mca_data[i]))!= 1) {
			fprintf(stderr,"Error reading DATA segment of %s at position %i\n",spe_file,i);
			return -1;
		}
	}

	fclose(filePtr);
	return 0;
}

double calculate_detector_absorption(struct xmi_input *input, int Z, int line) {
        double energy = LineEnergy(Z, line);
        int i, j;
        double rv = 1.0;

        for (i = 0 ; i < input->absorbers->n_det_layers ; i++) {
                double mu = 0.0;
                for (j = 0 ; j < input->absorbers->det_layers[i].n_elements ; j++)
			mu += CS_Total_Kissel(input->absorbers->det_layers[i].Z[j], energy) * input->absorbers->det_layers[i].weight[j];
                rv *= exp(-1.0*input->absorbers->det_layers[i].density * input->absorbers->det_layers[i].thickness * mu);
        }
        for (i = 0 ; i < input->detector->n_crystal_layers ; i++) {
                double mu = 0.0;
                for (j = 0 ; j < input->detector->crystal_layers[i].n_elements ; j++)
                        mu += CS_Total_Kissel(input->detector->crystal_layers[i].Z[j], energy) * input->detector->crystal_layers[i].weight[j];
                rv *= -1.0 * expm1(-1.0*input->detector->crystal_layers[i].density*input->detector->crystal_layers[i].thickness * mu);
        }
        return rv;
}


int main(int argc, char *argv[])
{
	char *xmim_input, *spe_file;
//	char *xmim_command;
//	char str[60];
//	FILE *fp;


	// argc expects 2 arguments: xmim input file and spe file
	if(argc != 3){
		printf("USAGE: 2 arguments should be supplied: XMSI and SPE file.\n");
		exit(0);
	}

	//First argument is XMSI file
	xmim_input = strdup(argv[1]);
	//Second argument is SPE file
	spe_file = strdup(argv[2]);

#if DEBUG == 1
	printf("XMSI: %s, SPE: %s.\n",xmim_input, spe_file);
#endif

	//HERE THE CODE OF xmim pymca.c BEGINS
	char *xmimsim_hdf5_solid_angles=NULL;
        struct xmi_input *xi = NULL;
//        struct xmi_pymca *xp = NULL ;
//        struct xmi_layer *matrix;
        double *weights_arr_quant;
        int i,j;
        xmi_inputFPtr inputFPtr;
        xmi_hdf5FPtr hdf5FPtr;
        GError *error = NULL;
        GOptionContext *context;
        gchar *hdf5_file=NULL;
	
	struct xmi_main_options xmi_options = xmi_get_default_main_options();
	struct xmi_main_options *options = &xmi_options;
        gchar *spe_file_noconv=NULL;
        gchar *spe_file_conv=NULL;
        gchar *csv_file_noconv=NULL;
        gchar *csv_file_conv=NULL;
        gchar *svg_file_noconv=NULL;
        gchar *svg_file_conv=NULL;
        gchar *excitation_file = NULL;
        gchar *geometry_file = NULL;
        gchar *detector_file = NULL;

	FILE *outPtr, *csv_convPtr, *csv_noconvPtr;
        gchar *filename;
        int use_rayleigh_normalization = 0;
        int use_roi_normalization = 0;
        int use_matrix_override= 0;
        int use_single_run= 0;
        int rayleigh_channel;
        double max_scale;
        double *scale, sum_scale, *k_exp, *k_sim, *l_exp, *l_sim;
        struct xmi_escape_ratios *escape_ratios_def=NULL;
        char *xmimsim_hdf5_escape_ratios = NULL;
//        double sum_roi;
        static int version = 0;
//        gchar *xmimsim_hdf5_solid_angles_utf8 = NULL;
//        gchar *xmimsim_hdf5_escape_ratios_utf8 = NULL;
//        gchar *hdf5_file_utf8 = NULL;

        GArray *entries = g_array_sized_new(TRUE, FALSE, sizeof(GOptionEntry), 30);
#define ADD_OPTION(long_name, short_name, flags, arg, arg_data, description, arg_description) \
        { \
                GOptionEntry entry = {long_name, short_name, flags, arg, arg_data, description, arg_description}; \
                g_array_append_val(entries, entry); \
        }
        ADD_OPTION("enable-M-lines", 0, 0, G_OPTION_ARG_NONE, &options->use_M_lines, "Enable M lines (default)", NULL );
        ADD_OPTION("disable-M-lines", 0, G_OPTION_FLAG_REVERSE, G_OPTION_ARG_NONE, &options->use_M_lines, "Disable M lines", NULL );
        ADD_OPTION("spe-file-unconvoluted",0,0,G_OPTION_ARG_FILENAME,&spe_file_noconv,"Write detector unconvoluted spectra to file",NULL);
        ADD_OPTION("spe-file",0,0,G_OPTION_ARG_FILENAME,&spe_file_conv,"Write detector convoluted spectra to file",NULL);
        ADD_OPTION("csv-file-unconvoluted",0,0,G_OPTION_ARG_FILENAME,&csv_file_noconv,"Write detector unconvoluted spectra to CSV file",NULL);
        ADD_OPTION("csv-file",0,0,G_OPTION_ARG_FILENAME,&csv_file_conv,"Write detector convoluted spectra to CSV file",NULL);
        ADD_OPTION("svg-file-unconvoluted",0,0,G_OPTION_ARG_FILENAME,&svg_file_noconv,"Write detector unconvoluted spectra to SVG file",NULL);
        ADD_OPTION("svg-file",0,0,G_OPTION_ARG_FILENAME,&svg_file_conv,"Write detector convoluted spectra to SVG file",NULL);
        ADD_OPTION("with-hdf5-data",0,0,G_OPTION_ARG_FILENAME,&hdf5_file,"Select a HDF5 data file (advanced usage)",NULL);
        ADD_OPTION("enable-scatter-normalization", 0, 0, G_OPTION_ARG_NONE,&use_rayleigh_normalization,"Enable Rayleigh peak based intensity normalization",NULL);
        ADD_OPTION("disable-scatter-normalization", 0, G_OPTION_FLAG_REVERSE, G_OPTION_ARG_NONE,&use_rayleigh_normalization,"Disable Rayleigh peak based intensity normalization (default)",NULL);
        ADD_OPTION("enable-roi-normalization", 0, 0, G_OPTION_ARG_NONE,&use_roi_normalization,"Enable region of interest integration based intensity normalization",NULL);
        ADD_OPTION("disable-roi-normalization", 0, G_OPTION_FLAG_REVERSE, G_OPTION_ARG_NONE,&use_roi_normalization,"Disable region of interest integration based intensity normalization (default)",NULL);
        ADD_OPTION("enable-matrix-override", 0, 0, G_OPTION_ARG_NONE,&use_matrix_override,"If the matrix includes quantifiable elements, use a similar matrix instead",NULL);
        ADD_OPTION("disable-matrix-override", 0, G_OPTION_FLAG_REVERSE, G_OPTION_ARG_NONE,&use_matrix_override,"If the matrix includes quantifiable elements, do not use a similar matrix instead (default)",NULL);
        ADD_OPTION("enable-pile-up", 0, 0, G_OPTION_ARG_NONE, &options->use_sum_peaks, "Enable pile-up", NULL );
        ADD_OPTION("disable-pile-up", 0, G_OPTION_FLAG_REVERSE, G_OPTION_ARG_NONE, &options->use_sum_peaks, "Disable pile-up (default)", NULL );
        ADD_OPTION("enable-escape-peaks", 0, 0, G_OPTION_ARG_NONE, &options->use_escape_peaks, "Enable escape peaks (default)", NULL );
        ADD_OPTION("disable-escape-peaks", 0, G_OPTION_FLAG_REVERSE, G_OPTION_ARG_NONE, &options->use_escape_peaks, "Disable escape peaks", NULL );
        ADD_OPTION("enable-poisson", 0, 0, G_OPTION_ARG_NONE, &options->use_poisson, "Generate Poisson noise in the spectra", NULL );
        ADD_OPTION("disable-poisson", 0, G_OPTION_FLAG_REVERSE, G_OPTION_ARG_NONE, &options->use_poisson, "Disable the generating of spectral Poisson noise (default)", NULL );
#if defined(HAVE_OPENCL_CL_H) || defined(HAVE_CL_CL_H)
        ADD_OPTION("enable-opencl", 0, 0, G_OPTION_ARG_NONE, &options->use_opencl, "Enable OpenCL (default)", NULL );
        ADD_OPTION("disable-opencl", 0, G_OPTION_FLAG_REVERSE, G_OPTION_ARG_NONE, &options->use_opencl, "Disable OpenCL", NULL );
#endif
        ADD_OPTION("enable-advanced-compton", 0, 0, G_OPTION_ARG_NONE, &options->use_advanced_compton, "Enable advanced yet slower Compton simulation", NULL );
        ADD_OPTION("disable-advanced-compton", 0, G_OPTION_FLAG_REVERSE, G_OPTION_ARG_NONE, &options->use_advanced_compton, "Disable advanced yet slower Compton simulation (default)", NULL );
        ADD_OPTION("enable-single-run", 0, 0, G_OPTION_ARG_NONE, &use_single_run, "Force the simulation to run just once", NULL );
        ADD_OPTION("override-excitation",0,0,G_OPTION_ARG_FILENAME,&excitation_file, "Override excitation from XMSI file",NULL);
        ADD_OPTION("override-detector",0,0,G_OPTION_ARG_FILENAME,&detector_file, "Override detector from XMSI file",NULL);
        ADD_OPTION("override-geometry",0,0,G_OPTION_ARG_FILENAME,&geometry_file, "Override geometry from XMSI file",NULL);
        ADD_OPTION("set-threads",0,0,G_OPTION_ARG_INT, &options->omp_num_threads, "Sets the number of threads to NTHREADS (default=max)", "NTHREADS");
        ADD_OPTION("verbose", 'v', 0, G_OPTION_ARG_NONE, &options->verbose, "Verbose mode", NULL );
        ADD_OPTION("version", 0, 0, G_OPTION_ARG_NONE, &version, "Display version information", NULL );

//#if DEBUG == 1
	options->verbose = 1;
	options->use_sum_peaks = 1;
//#endif

	double *channels;
        double **channels_conv;
        double zero_sum;
        double *brute_history;
        double *var_red_history;
        double sum_k, sum_l, sum_temp, sum_weights;
        struct xmi_solid_angle *solid_angle_def;
        char *xmi_input_string;


        xmi_init_hdf5();
        setbuf(stdout,NULL);
        //locale...
        //        //setlocale(LC_ALL,"C");
        g_setenv("LANG","en_US",TRUE);

        //parse options
        context = g_option_context_new ("inputfile outputfile");
        g_option_context_add_main_entries(context, (const GOptionEntry *) entries->data, NULL);
        g_option_context_set_summary(context, "xmimsim-pymca: a program for the quantification of X-ray fluorescence spectra using inverse Monte-Carlo simulations. Inputfiles should be prepared using PyMCA\n");
        if (!g_option_context_parse (context, &argc, &argv, &error)) {
                g_fprintf(stdout, "option parsing failed: %s\n", error->message);
                return 1;
        }

        if (version) {
                g_fprintf(stdout,"%s",xmi_version_string());
                return 0;
        }

	//check for conflicting options
	if (use_rayleigh_normalization + use_roi_normalization + use_matrix_override + use_single_run > 1) {
                g_fprintf(stdout,"Options conflict: Use either --enable-rayleigh-normalization or --enable-roi-normalization or --enable-matrix-override or --enable-single-run. No combinations of these are allowed\n");
                exit(1);
        }

        if (options->omp_num_threads > xmi_omp_get_max_threads() ||
                        options->omp_num_threads < 1) {
                options->omp_num_threads = xmi_omp_get_max_threads();
        }

	//load xml catalog
	if (xmi_xmlLoadCatalog() == 0) {
                g_fprintf(stdout, "Could not load catalog: %s\n", error->message);
                return 1;
        }
        else if (options->verbose)
                g_fprintf(stdout,"XML catalog loaded\n");

	if (xmi_get_hdf5_data_file(&hdf5_file) == 0){
		return 1;
	}

	//start random number acquisition
	if (xmi_start_random_acquisition() == 0){
		return 1;
	}


	//allocate and generate the default input file, then read in actual input file
        xmi_read_input_xml(xmim_input, &xi);
#if DEBUG == 1
	xmi_print_input(stdout,xi);
#endif

	//Read in SPE data
	struct xmo_spe_data *spe;

	spe = malloc(sizeof(struct xmo_spe_data));
	if (xmo_read_spe(spe_file, spe) != 0) {
		g_fprintf(stdout,"Could not read spe file %s\n",spe_file);
	}
#if DEBUG == 1
	printf("SPE:%s NCHAN: %i\n",spe_file, spe->nchannels);
#endif




	//calculate initial
	//	Just have program use elements already in xmsi input file... User should adjust this
#define XMO_N_Z_QUANT xi->composition->layers[xi->composition->reference_layer-1].n_elements
#define XMO_SCAT_INT (xi->excitation->discrete[xi->excitation->n_discrete-1].horizontal_intensity+xi->excitation->discrete[xi->excitation->n_discrete-1].vertical_intensity)
	weights_arr_quant = (double *) g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	for (i = 0; i < XMO_N_Z_QUANT; i++){
		weights_arr_quant[i] = xi->composition->layers[xi->composition->reference_layer-1].weight[i];
#if DEBUG == 1
printf("**%i, %lf\n", xi->composition->layers[xi->composition->reference_layer-1].Z[i], xi->composition->layers[xi->composition->reference_layer-1].weight[i] );
#endif
	}

	//copy to the corresponding fortran variable
	xmi_input_C2F(xi,&inputFPtr);
	//initialization
	if (xmi_init_input(&inputFPtr) == 0) {
		return 1;
	}
	//initialize HDF5 data
	//read from HDF5 file what needs to be read in
	if (xmi_init_from_hdf5(hdf5_file,inputFPtr,&hdf5FPtr,*options) == 0) {
		g_fprintf(stdout,"Could not initialize from hdf5 data file\n");
		return 1;
	}
	else g_fprintf(stdout,"HDF5 datafile %s successfully processed\n",hdf5_file);

	xmi_update_input_from_hdf5(inputFPtr, hdf5FPtr);

	//determine filename first
	if (xmi_get_solid_angle_file(&xmimsim_hdf5_solid_angles, 1) == 0)
		return 1;
	g_fprintf(stdout,"Querying %s for solid angle grid\n",xmimsim_hdf5_solid_angles);

	
	//check if solid angles are already precalculated
	if (xmi_find_solid_angle_match(xmimsim_hdf5_solid_angles , xi, &solid_angle_def, *options) == 0){
#if DEBUG == 1
		printf("xmi_find_solid_angle_match == 0\n");
#endif
		return 1;
	}
	if (solid_angle_def == NULL) {
		g_fprintf(stdout,"Precalculating solid angle grid\n");
		//doesn't exist yet
		//convert input to string
		if (xmi_write_input_xml_to_string(&xmi_input_string, xi) == 0) {
			g_fprintf(stdout, "Could not write input to XML string: %s\n",xmi_input_string);
			return 1;
		}
		xmi_solid_angle_calculation(inputFPtr, &solid_angle_def, xmi_input_string, *options);
		g_free(xmi_input_string);
		//update hdf5 file
		if( xmi_update_solid_angle_hdf5_file(xmimsim_hdf5_solid_angles , solid_angle_def) == 0)
			return 1;
		else g_fprintf(stdout,"%s was successfully updated with new solid angle grid\n",xmimsim_hdf5_solid_angles);
	}
	else g_fprintf(stdout,"Solid angle grid already present in %s\n",xmimsim_hdf5_solid_angles);

	//everything is read in... start iteration
	i = 0;
	//convolute_spectrum
	channels_conv = (double **) g_malloc0(sizeof(double *)*(xi->general->n_interactions_trajectory+1));

#define ARRAY3D_FORTRAN(array,i,j,k,Ni,Nj,Nk) (array[Nj*Nk*(i-1)+Nk*(j-1)+(k-1)])

	sum_k = sum_l = XMO_CONV_THRESHOLD*10.0;

	channels = NULL;
	brute_history = NULL;
	var_red_history = NULL;

	//read escape ratios 
	if (options->use_escape_peaks) {
		if (xmi_get_escape_ratios_file(&xmimsim_hdf5_escape_ratios, 1) == 0)
			return 1;

		if (options->verbose)
			g_fprintf(stdout,"Querying %s for escape peak ratios\n",xmimsim_hdf5_escape_ratios);


		//check if escape ratios are already precalculated
		if (xmi_find_escape_ratios_match(xmimsim_hdf5_escape_ratios , xi, &escape_ratios_def, *options) == 0)
			return 1;
		if (escape_ratios_def == NULL) {
			if (options->verbose)
				g_fprintf(stdout,"Precalculating escape peak ratios\n");
			//doesn't exist yet
			//convert input to string
			if (xmi_write_input_xml_to_string(&xmi_input_string, xi) == 0) {
				g_fprintf(stdout, "Could not write input to XML string: %s\n", error->message);
				return 1;
			}
			xmi_escape_ratios_calculation(xi, &escape_ratios_def, xmi_input_string,hdf5_file,*options, xmi_get_default_escape_ratios_options());
			//update hdf5 file
			if( xmi_update_escape_ratios_hdf5_file(xmimsim_hdf5_escape_ratios , escape_ratios_def) == 0)
				return 1;
			else if (options->verbose)
				g_fprintf(stdout,"%s was successfully updated with new escape peak ratios\n",xmimsim_hdf5_escape_ratios);
		}
		else if (options->verbose)
			g_fprintf(stdout,"Escape peak ratios already present in %s\n",xmimsim_hdf5_escape_ratios);
	}
	else if (options->verbose)
		g_fprintf(stdout,"No escape peaks requested: calculation is redundant\n");

	if (use_rayleigh_normalization && xi->excitation->discrete[xi->excitation->n_discrete-1].energy > 0.0 &&  XMO_SCAT_INT > 0.0) {
		rayleigh_channel = (int) ((xi->excitation->discrete[xi->excitation->n_discrete-1].energy - xi->detector->zero)/xi->detector->gain);
		if (rayleigh_channel > xi->detector->nchannels) {
			g_fprintf(stdout,"Channel of excitation energy is not included in spectrum from pymca\n");
			return 1;
		}
#if DEBUG == 1
		g_fprintf(stdout,"rayleigh_channel: %i\n", rayleigh_channel);
#endif
	}
	
	scale = (double *) g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	k_exp= (double *) g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	k_sim = (double *) g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	l_exp = (double *) g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	l_sim= (double *) g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	

	//fill up k_exp and l_exp
	//contains the intensities of the experimental spectrum
	double *energy_k, *energy_l; 
	double *k_sim_low, *l_sim_low, *k_sim_high, *l_sim_high, *weight_low, *weight_high;
	energy_k = g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	energy_l = g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	k_sim_low = g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	l_sim_low = g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	k_sim_high = g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	l_sim_high = g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	weight_low = g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	weight_high = g_malloc(sizeof(double)*XMO_N_Z_QUANT);
	for (j = 0 ; j < XMO_N_Z_QUANT ; j++) {
		energy_k[j] = LineEnergy(xi->composition->layers[xi->composition->reference_layer-1].Z[j], KA1_LINE);
		energy_l[j] = LineEnergy(xi->composition->layers[xi->composition->reference_layer-1].Z[j], LA1_LINE);
		k_sim_low[j] = 0.0;
		l_sim_low[j] = 0.0;
		k_sim_high[j] = DBL_MAX;
		l_sim_high[j] = DBL_MAX;
		weight_low[j] = 0.0;
		weight_high[j] = 1.1; //110% w%... simply an impossibly high concentration
#if DEBUG == 1
	printf("--Energy KA1: %lf LA1: %lf keV\n",energy_k[j],energy_l[j]);
#endif
		if (energy_k[j] > 1. && (energy_k[j] - xi->detector->zero)/xi->detector->gain < xi->detector->nchannels){
			k_exp[j] = (double)spe->mca_data[(int)((energy_k[j]-xi->detector->zero)/xi->detector->gain)];
			l_exp[j] = 0.0;
		}
		else if (energy_l[j] > 1. && (energy_l[j] - xi->detector->zero)/xi->detector->gain < xi->detector->nchannels){
			k_exp[j] = 0.0;
			l_exp[j] = (double)spe->mca_data[(int)((energy_l[j]-xi->detector->zero)/xi->detector->gain)];
		}
		else {
			k_exp[j] = 0.0;
			l_exp[j] = 0.0;
		}
	}

//	from xmi_pymca.c, perhaps not necessary here
//	if (use_roi_normalization) {
//		sum_roi = 5*xp->sum_xmin_xmax;
//	}
//	else {
//		sum_roi = xp->sum_xmin_xmax;
//	}

#define ARRAY2D_FORTRAN(array,i,j,Ni,Nj) (array[(Nj)*(i)+(j)])
	
	while ((sum_k > XMO_CONV_THRESHOLD) || (sum_l > XMO_CONV_THRESHOLD)){ //|| fabs(sum_roi-xp->sum_xmin_xmax)/xp->sum_xmin_xmax > 0.05
		xmi_deallocate(channels);
		xmi_deallocate(brute_history);
		xmi_deallocate(var_red_history);

		if (i++ > XMI_PYMCA_MAX_ITERATIONS) {
			g_fprintf(stdout,"No convergence after %i iterations... Fatal error\n",i);
			return 0;
		}
#if DEBUG == 1
		sprintf(tempFile, "xmimsim-pymca_debug_%i.xmsi",i);
		xmi_write_input_xml(tempFile, xi, NULL);
#endif



		//launch simulation
		if (options->verbose)
			g_fprintf(stdout,"Simulating interactions\n");
//if (i == 1) {
		if (xmi_main_msim(inputFPtr, hdf5FPtr, 1, &channels, *options, &brute_history, &var_red_history, solid_angle_def) == 0) {
			g_fprintf(stdout,"Error in xmi_main_msim\n");
			return 1;
		}
//}
		if (options->verbose)
			g_fprintf(stdout,"Interactions simulation finished\n");
#if DEBUG == 1
		//write input structure
		//xmi_print_input(stdout,xi);
		zero_sum = xmi_sum_double(channels, xi->detector->nchannels);
		//convolute_spectrum
        	double *channels_conv_temp2;
		channels_conv_temp = (double **) g_malloc0(sizeof(double *)*(xi->general->n_interactions_trajectory+1));

		for (j=(zero_sum > 0.0 ? 0 : 1) ; j <= xi->general->n_interactions_trajectory ; j++) {
			xmi_detector_convolute(inputFPtr, hdf5FPtr, channels+j*xi->detector->nchannels, &channels_conv_temp2, xi->detector->nchannels, options);
			channels_conv_temp[i] = xmi_memdup(channels_conv_temp2,sizeof(double)*xi->detector->nchannels);
		}
		//write to xml outputfile
		sprintf(tempFile, "xmimsim-pymca_debug_%i.xmso",i);
		if (xmi_write_output_xml(tempFile, xi, brute_history, options->use_variance_reduction == 1 ? var_red_history : NULL, channels_conv_temp, channels, xi->detector->nchannels, xmim_input, zero_sum > 0.0 ? 1 : 0, NULL) == 0) {
			return 1;
		}
		free(channels_conv_temp);
		free(channels_conv_temp2);
#endif

		//optimize concentrations
		//if normalization is enabled -> do not optimize after first run. Only the intensity of the exciting radiation will be adjusted in this case
		if (!(use_rayleigh_normalization && xi->excitation->discrete[xi->excitation->n_discrete-1].energy > 0.0 &&  XMO_SCAT_INT > 0.0 && i == 1) && !(use_roi_normalization && i % 2 == 1)) {
			if (options->verbose)
				g_fprintf(stdout, "Recalculating weight fractions\n");

			sum_k = sum_l = 0.0;
			sum_scale = 0.0;
			sum_weights = 0.0;

			//convolute simulated spectrum so we can obtain channel intensities
			double **channels_def_ptrs = g_malloc0(sizeof(double *) * (xi->general->n_interactions_trajectory+1));
			for (j = 0 ; j <= xi->general->n_interactions_trajectory ; j++)
				channels_def_ptrs[j] = channels+j*xi->detector->nchannels;
			xmi_detector_convolute_all(inputFPtr, channels_def_ptrs, channels_conv, brute_history, var_red_history, *options, escape_ratios_def, xi->general->n_interactions_trajectory, 1);
			g_free(channels_def_ptrs);


			for (j = XMO_N_Z_QUANT-1 ; j > -1 ; j--) { //TODO: it would be nicest to do this iteration from highest energy to lowest?

				k_sim[j] = 0.0;
				l_sim[j] = 0.0;

				if (energy_k[j] > 1. && (energy_k[j] - xi->detector->zero)/xi->detector->gain < xi->detector->nchannels){
					k_sim[j] = channels_conv[xi->general->n_interactions_trajectory][(int)((energy_k[j]-xi->detector->zero)/xi->detector->gain)];
					l_sim[j] = 0.0;
				}
				else if (energy_l[j] > 1. && (energy_l[j] - xi->detector->zero)/xi->detector->gain < xi->detector->nchannels){
					k_sim[j] = 0.0;
					l_sim[j] = channels_conv[xi->general->n_interactions_trajectory][(int)((energy_l[j]-xi->detector->zero)/xi->detector->gain)];
				}
				else {
					k_sim[j] = 0.0;
					l_sim[j] = 0.0;
				}

				// TODO: apply detector absorption to var_red values
//				for (k = 0 ; k <= xi->general->n_interactions_trajectory ; k++) {
//					k_sim[j] += ARRAY3D_FORTRAN(var_red_history, xi->composition->layers[xi->composition->reference_layer-1].Z[j], abs(KL2_LINE),k,100,385,xi->general->n_interactions_trajectory) * calculate_detector_absorption(xi, xi->composition->layers[xi->composition->reference_layer-1].Z[j], KL2_LINE);
//					k_sim[j] += ARRAY3D_FORTRAN(var_red_history, xi->composition->layers[xi->composition->reference_layer-1].Z[j], abs(KL3_LINE),k,100,385,xi->general->n_interactions_trajectory) * calculate_detector_absorption(xi, xi->composition->layers[xi->composition->reference_layer-1].Z[j], KL3_LINE);
//				}
//				for (k = 0 ; k <= xi->general->n_interactions_trajectory ; k++) {
//					l_sim[j] += ARRAY3D_FORTRAN(var_red_history, xi->composition->layers[xi->composition->reference_layer-1].Z[j], abs(L3M4_LINE),k,100,385,xi->general->n_interactions_trajectory) * calculate_detector_absorption(xi, xi->composition->layers[xi->composition->reference_layer-1].Z[j], L3M4_LINE);
//					l_sim[j] += ARRAY3D_FORTRAN(var_red_history, xi->composition->layers[xi->composition->reference_layer-1].Z[j], abs(L3M5_LINE),k,100,385,xi->general->n_interactions_trajectory) * calculate_detector_absorption(xi, xi->composition->layers[xi->composition->reference_layer-1].Z[j], L3M5_LINE);
//				}

//				if (k_exp[j] > 1.0 && k_sim[j] > 0.0) {
//					scale[j] = k_exp[j]/k_sim[j];
//					if (fabs(scale[j]-1.0) < 0.05 ) scale[j] = 1.0;
//				}
//				else if (l_exp[j] > 1.0 && l_sim[j] > 0.0) {
//					scale[j] = l_exp[j]/l_sim[j];
//					if (fabs(scale[j]-1.0) < 0.05 ) scale[j] = 1.0;
//				}
//				else scale[j] = 1.0;
//				sum_scale += scale[j]*weights_arr_quant[j];

				//do not allow too large jumps!
				//if weight <= 0.0001 then max is 100
				//else if weight <= 0.01 then max is 10
				//else if weight <= 0.1 then max is 2.5
				//else if weight <= 0.25 then max is 1.5
				//else if weight <= 0.375 then max is 1.20
				if (weights_arr_quant[j] <= 0.0001)
					max_scale = 100.0;
				else if (weights_arr_quant[j] <= 0.01)
					max_scale = 10.0;
				else if (weights_arr_quant[j] <= 0.1)
					max_scale = 2.5;
				else if (weights_arr_quant[j] <= 0.25)
					max_scale = 1.5;
				else if (weights_arr_quant[j] <= 0.375)
					max_scale = 1.2;
				else if (weights_arr_quant[j] <= 0.50)
					max_scale = 1.1;
				else if (weights_arr_quant[j] <= 0.6)
					max_scale = 1.05;
				else if (weights_arr_quant[j] <= 0.7)
					max_scale = 1.025;
				else
					max_scale = 1.01;

				// store simulated intensity and corresponding concentrations
				// and calculate scaling factor
				if (k_sim[j] < k_exp[j] && k_sim[j] > 0.0 && k_exp[j] > 1.0) {
#if DEBUG == 1
	fprintf(stdout,"**Case1");
#endif
					if (k_sim[j] > k_sim_low[j]) {
#if DEBUG == 1
fprintf(stdout,"-");
#endif
						k_sim_low[j] = k_sim[j];
						weight_low[j] = weights_arr_quant[j];
					}
					// interpolate between known low and high value to obtain new scale estimate
					if (weight_low[j] > 0.0 && weight_high[j] < 1.) {
						scale[j] = ( (k_exp[j]-k_sim_low[j]) / ((k_sim_high[j]-k_sim_low[j])/(weight_high[j]-weight_low[j])) + weight_low[j]) / weights_arr_quant[j];
						//we should converge for this peak, yet it could be influenced by secondary excitation etc
						//	so reset low and high values
						k_sim_low[j] = 0.0;
						k_sim_high[j] = DBL_MAX;
						weight_low[j] = 0.0;
						weight_high[j] = 1.1;
#if DEBUG == 1
	fprintf(stdout,"a\n");
#endif
					}
					else {
						scale[j] = k_exp[j]/k_sim[j];
#if DEBUG == 1
	fprintf(stdout,"b\n");
#endif
					}
				}
				else if (l_sim[j] < l_exp[j] && l_sim[j] > 0.0 && l_exp[j] > 1.0) {
#if DEBUG == 1
	fprintf(stdout,"**Case2");
#endif
					if (l_sim[j] > l_sim_low[j]) {
#if DEBUG == 1
	fprintf(stdout,"-");
#endif
						l_sim_low[j] = l_sim[j];
						weight_low[j] = weights_arr_quant[j];
					}
					// interpolate between known low and high value to obtain new scale estimate
					if (weight_low[j] > 0.0 && weight_high[j] < 1.) {
						scale[j] = ( (l_exp[j]-l_sim_low[j]) / ((l_sim_high[j]-l_sim_low[j])/(weight_high[j]-weight_low[j])) + weight_low[j]) / weights_arr_quant[j];
						//we should converge for this peak, yet it could be influenced by secondary excitation etc
						//	so reset low and high values
						l_sim_low[j] = 0.0;
						l_sim_high[j] = DBL_MAX;
						weight_low[j] = 0.0;
						weight_high[j] = 1.1;
#if DEBUG == 1
	fprintf(stdout,"a\n");
#endif
					}
					else {
						scale[j] = l_exp[j]/l_sim[j];
#if DEBUG == 1
	fprintf(stdout,"b\n");
#endif
					}
				}
				else if (k_sim[j] > k_exp[j] && k_sim[j] > 0.0 && k_exp[j] > 1.0) {
#if DEBUG == 1
	fprintf(stdout,"**Case3");
#endif
					if (k_sim[j] < k_sim_high[j]) {
#if DEBUG == 1
	fprintf(stdout,"-");
#endif
						k_sim_high[j] = k_sim[j];
						weight_high[j] = weights_arr_quant[j];
					}
					// interpolate between known low and high value to obtain new scale estimate
					if (weight_low[j] > 0.0 && weight_high[j] < 1.) {
						scale[j] = ( (k_exp[j]-k_sim_low[j]) / ((k_sim_high[j]-k_sim_low[j])/(weight_high[j]-weight_low[j])) + weight_low[j]) / weights_arr_quant[j];
						//we should converge for this peak, yet it could be influenced by secondary excitation etc
						//	so reset low and high values
						k_sim_low[j] = 0.0;
						k_sim_high[j] = DBL_MAX;
						weight_low[j] = 0.0;
						weight_high[j] = 1.1;
#if DEBUG == 1
	fprintf(stdout,"a\n");
#endif
					}
					else {
						scale[j] = k_exp[j]/k_sim[j];
#if DEBUG == 1
	fprintf(stdout,"b\n");
#endif
					}
				}
				else if (l_sim[j] > l_exp[j] && l_sim[j] > 0.0 && l_exp[j] > 1.0) {
#if DEBUG == 1
	fprintf(stdout,"**Case4");
#endif
					if (l_sim[j] < l_sim_high[j]) {
#if DEBUG == 1
	fprintf(stdout,"-");
#endif
						l_sim_high[j] = l_sim[j];
						weight_high[j] = weights_arr_quant[j];
					}
					// interpolate between known low and high value to obtain new scale estimate
					if (weight_low[j] > 0.0 && weight_high[j] < 1.) {
						scale[j] = ( (l_exp[j]-l_sim_low[j]) / ((l_sim_high[j]-l_sim_low[j])/(weight_high[j]-weight_low[j])) + weight_low[j]) / weights_arr_quant[j];
						//we should converge for this peak, yet it could be influenced by secondary excitation etc
						//	so reset low and high values
						l_sim_low[j] = 0.0;
						l_sim_high[j] = DBL_MAX;
						weight_low[j] = 0.0;
						weight_high[j] = 1.1;
#if DEBUG == 1
	fprintf(stdout,"a\n");
#endif
					}
					else {
						scale[j] = l_exp[j]/l_sim[j];
#if DEBUG == 1
	fprintf(stdout,"b\n");
#endif
					}
				}
				else {
					scale[j] = 1.0;
#if DEBUG == 1
	fprintf(stdout,"**Case5\n");
#endif
				}

				if (fabs(k_sim[j]-k_exp[j]) < sqrt(k_sim[j]) && k_sim[j] > 0.0 && k_exp[j] > 1.0) {
					scale[j] = 1.0;
					//we have converged for this peak, yet it could be influenced by secondary excitation etc
					//	so reset low and high values
					k_sim_low[j] = 0.0;
					k_sim_high[j] = DBL_MAX;
					weight_low[j] = 0.0;
					weight_high[j] = 1.1;
#if DEBUG == 1
	fprintf(stdout,"**Case6\n");
#endif
				}
				else if (fabs(l_sim[j]-l_sim[j]) < sqrt(l_sim[j]) && l_sim[j] > 0.0 && l_exp[j] > 1.0) {
					scale[j] = 1.0;
					//we have converged for this peak, yet it could be influenced by secondary excitation etc
					//	so reset low and high values
					l_sim_low[j] = 0.0;
					l_sim_high[j] = DBL_MAX;
					weight_low[j] = 0.0;
					weight_high[j] = 1.1;
#if DEBUG == 1
	fprintf(stdout,"**Case7\n");
#endif
				}
				sum_scale += scale[j]*weights_arr_quant[j]; //TODO: probably safe to remove this line and variable sum_scale


				if (scale[j] > max_scale) {
					weights_arr_quant[j] *=	max_scale;
				}
				else {
					weights_arr_quant[j] *= scale[j];
				}
				if (k_exp[j] > 1.0 && k_sim[j] > 0.0  ) {
					sum_temp = 0.0;
					sum_temp += (k_exp[j]-k_sim[j])*(k_exp[j]-k_sim[j]);
					sum_temp /= k_exp[j];
					sum_temp /= k_exp[j];
					sum_k += sum_temp;
				}
				else if (l_exp[j] > 1.0 && l_sim[j] > 0.0  ) {
					sum_temp = 0.0;
					sum_temp += (l_exp[j]-l_sim[j])*(l_exp[j]-l_sim[j]);
					sum_temp /= l_exp[j];
					sum_temp /= l_exp[j];
					sum_l += sum_temp;
				}
				//K-lines
				//make history_sum

//#if DEBUG == 1
				fprintf(stdout,"Element :%i\n",xi->composition->layers[xi->composition->reference_layer-1].Z[j]);
				if (k_sim[j] != 0.0) {
					fprintf(stdout,"k_exp[j]: %lf\n",k_exp[j]);
					fprintf(stdout,"k_sim[j]: %lf\n",k_sim[j]);
				}
				else if (l_sim[j] != 0.0) {
					fprintf(stdout,"l_exp[j]: %lf\n",l_exp[j]);
					fprintf(stdout,"l_sim[j]: %lf\n",l_sim[j]);
				}
				fprintf(stdout,"scale[j]: %lf\n",scale[j]);
				fprintf(stdout,"weight[j]: %lf\n",weights_arr_quant[j]);
	//			fprintf(stdout,"scatter_intensity from file: %lf\n",xp->scatter_intensity);
	//			fprintf(stdout,"scatter_intensity from MC: %lf\n",ARRAY2D_FORTRAN(channels,xi->general->n_interactions_trajectory,rayleigh_channel,xi->general->n_interactions_trajectory+1,xp->nchannels));
//#endif

			}

//			for (j = 0 ; j < XMO_N_Z_QUANT ; j++) {
//				//do not allow too large jumps!
//				//if weight <= 0.0001 then max is 100
//				//else if weight <= 0.01 then max is 10
//				//else if weight <= 0.1 then max is 2.5
//				//else if weight <= 0.25 then max is 1.5
//				//else if weight <= 0.375 then max is 1.20
//
//				//scale[j] /= sum_scale;
//
//				if (weights_arr_quant[j] <= 0.0001)
//					max_scale = 100.0;
//				else if (weights_arr_quant[j] <= 0.01)
//					max_scale = 10.0;
//				else if (weights_arr_quant[j] <= 0.1)
//					max_scale = 2.5;
//				else if (weights_arr_quant[j] <= 0.25)
//					max_scale = 1.5;
//				else if (weights_arr_quant[j] <= 0.375)
//					max_scale = 1.2;
//				else if (weights_arr_quant[j] <= 0.50)
//					max_scale = 1.1;
//				else if (weights_arr_quant[j] <= 0.6)
//					max_scale = 1.05;
//				else if (weights_arr_quant[j] <= 0.7)
//					max_scale = 1.025;
//				else
//					max_scale = 1.01;
//
//				if (k_exp[j] > 1.0 && k_sim[j] > 0.0  ) {
//					sum_temp = 0.0;
//					sum_temp += (k_exp[j]-k_sim[j])*(k_exp[j]-k_sim[j]);
//					sum_temp /= k_exp[j];
//					sum_temp /= k_exp[j];
//					sum_k += sum_temp;
//					if (scale[j] > max_scale) {
//						weights_arr_quant[j] *=	max_scale;
//					}
//					else {
//						weights_arr_quant[j] *= scale[j];
//					}
//				}
//				else if (l_exp[j] > 1.0 && l_sim[j] > 0.0  ) {
//					sum_temp = 0.0;
//					sum_temp += (l_exp[j]-l_sim[j])*(l_exp[j]-l_sim[j]);
//					sum_temp /= l_exp[j];
//					sum_temp /= l_exp[j];
//					sum_l += sum_temp;
//					if (scale[j] > max_scale) {
//						weights_arr_quant[j] *=	max_scale;
//					}
//					else {
//						weights_arr_quant[j] *= scale[j];
//					}
//				}
//#if DEBUG == 1
//			fprintf(stdout,"Element :%i\n",xi->composition->layers[xi->composition->reference_layer-1].Z[j]);
//			fprintf(stdout,"k_exp[j]: %lf\n",k_exp[j]);
//			fprintf(stdout,"k_sim[j]: %lf\n",k_sim[j]);
//			fprintf(stdout,"l_exp[j]: %lf\n",l_exp[j]);
//			fprintf(stdout,"l_sim[j]: %lf\n",l_sim[j]);
//			fprintf(stdout,"scale[j]: %lf\n",scale[j]);
//			fprintf(stdout,"weight[j]: %lf\n",weights_arr_quant[j]);
//#endif
//			}
			//Normalize weights
			for (j = 0 ; j < XMO_N_Z_QUANT ; j++) {
				sum_weights += weights_arr_quant[j];
			}
		}
		//update energy intensities when required
		if (use_rayleigh_normalization && xi->excitation->discrete[xi->excitation->n_discrete-1].energy > 0.0 && XMO_SCAT_INT > 0.0) {
			if (options->verbose)
				g_fprintf(stdout, "Scaling beam intensity according to Rayleigh signal\n");
			for (j = 0 ; j < xi->excitation->n_discrete ; j++) {
				xi->excitation->discrete[j].horizontal_intensity *= XMO_SCAT_INT/ARRAY2D_FORTRAN(channels,xi->general->n_interactions_trajectory,rayleigh_channel,xi->general->n_interactions_trajectory+1,xi->detector->nchannels);
				xi->excitation->discrete[j].vertical_intensity *= XMO_SCAT_INT/ARRAY2D_FORTRAN(channels,xi->general->n_interactions_trajectory,rayleigh_channel,xi->general->n_interactions_trajectory+1,xi->detector->nchannels);
			}

			for (j = 0 ; j < xi->excitation->n_continuous ; j++) {
				xi->excitation->continuous[j].horizontal_intensity *= XMO_SCAT_INT/ARRAY2D_FORTRAN(channels,xi->general->n_interactions_trajectory,rayleigh_channel,xi->general->n_interactions_trajectory+1,xi->detector->nchannels);
				xi->excitation->continuous[j].vertical_intensity *= XMO_SCAT_INT/ARRAY2D_FORTRAN(channels,xi->general->n_interactions_trajectory,rayleigh_channel,xi->general->n_interactions_trajectory+1,xi->detector->nchannels);
			}
			if (i > 1) {
				//update concentrations in input
				xi->composition->layers[xi->composition->reference_layer-1].weight[i] = weights_arr_quant[i]/sum_weights;
				//g_free(xi->composition->layers[xp->ilay_pymca].Z);
				//g_free(xi->composition->layers[xp->ilay_pymca].weight);
				//xi->composition->layers[xp->ilay_pymca] = xmi_ilay_composition_pymca(matrix, xp, weights_arr_quant);
			}
		}
		else if (use_roi_normalization) {
			//integrate the region of interest
			//first apply detector convolution
			if (i % 2 == 1) {
				if (options->verbose)
					g_fprintf(stdout, "Scaling beam intensity according to region of interest intensity integration\n");
//        			double *channels_conv_temp2;
//				xmi_detector_convolute_spectrum(inputFPtr, channels+xi->general->n_interactions_trajectory*xi->detector->nchannels, &channels_conv_temp2, options, escape_ratios_def, xi->general->n_interactions_trajectory);

//				sum_roi = 0.0;
//				for (j = xp->xmin ; j <= xp->xmax ; j++)
//					sum_roi += channels_conv_temp2[j];

//				g_fprintf(stdout,"sum_roi: %lf\n", sum_roi);

//				xmi_deallocate(channels_conv_temp2);

//				for (j = 0 ; j < xi->excitation->n_discrete ; j++) {
//					xi->excitation->discrete[j].horizontal_intensity *= xp->sum_xmin_xmax/sum_roi;
//					xi->excitation->discrete[j].vertical_intensity *= xp->sum_xmin_xmax/sum_roi;
//				}

//				for (j = 0 ; j < xi->excitation->n_continuous ; j++) {
//					xi->excitation->continuous[j].horizontal_intensity *= xp->sum_xmin_xmax/sum_roi;
//					xi->excitation->continuous[j].vertical_intensity *= xp->sum_xmin_xmax/sum_roi;
//				}
			}
			else if (i % 2 == 0) {
				//update concentrations in input
				xi->composition->layers[xi->composition->reference_layer-1].weight[i] = weights_arr_quant[i]/sum_weights;
				//g_free(xi->composition->layers[xp->ilay_pymca].Z);
				//g_free(xi->composition->layers[xp->ilay_pymca].weight);
				//xi->composition->layers[xp->ilay_pymca] = xmi_ilay_composition_pymca(matrix, xp, weights_arr_quant);
			}
		}
		else {
			//update concentrations in input
			xi->composition->layers[xi->composition->reference_layer-1].weight[i] = weights_arr_quant[i]/sum_weights;
			//g_free(xi->composition->layers[xp->ilay_pymca].Z);
			//g_free(xi->composition->layers[xp->ilay_pymca].weight);
			//xi->composition->layers[xp->ilay_pymca] = xmi_ilay_composition_pymca(matrix, xp, weights_arr_quant);
		}

		//reload fortran input
		xmi_free_input_F(&inputFPtr);
		xmi_input_C2F(xi, &inputFPtr);
		if (xmi_init_input(&inputFPtr) == 0) {
			return 1;
		}
		xmi_update_input_from_hdf5(inputFPtr, hdf5FPtr);

		g_fprintf(stdout,"Iteration: %i\n",i);
		g_fprintf(stdout,"sum_k: %g\n",sum_k);
		g_fprintf(stdout,"sum_l: %g\n",sum_l);
	} //iterative while loop

	xmi_free_hdf5_F(&hdf5FPtr);

	zero_sum = xmi_sum_double(channels, xi->detector->nchannels);
//	//convolute_spectrum
//	channels_conv = (double **) g_malloc0(sizeof(double *)*(xi->general->n_interactions_trajectory+1));

	double **channels_def_ptrs = g_malloc0(sizeof(double *) * (xi->general->n_interactions_trajectory+1));
	for (i = 0 ; i <= xi->general->n_interactions_trajectory ; i++)
		channels_def_ptrs[i] = channels+i*xi->detector->nchannels;

	
	xmi_detector_convolute_all(inputFPtr, channels_def_ptrs, channels_conv, brute_history, var_red_history, *options, escape_ratios_def, xi->general->n_interactions_trajectory, 1);

	g_free(channels_def_ptrs);

	if (xmi_end_random_acquisition() == 0) {
		return 1;
	}

	//write to xml outputfile
	struct xmi_output *output = xmi_output_raw2struct(xi, brute_history, var_red_history, channels_conv, channels, xmim_input, 1);
	if (xmi_write_output_xml(argv[2], output) == 0) {
		g_fprintf(stdout, "Could not write to %s: %s\n", argv[2], error->message);
		return 1;
	}
	else g_fprintf(stdout,"Output written to XMSO file %s\n",argv[argc-1]); //TODO: should make sure this is the input file ending with .xmso
	xmi_free_output(output);

	//write to CSV and SPE if necessary...
	csv_convPtr = csv_noconvPtr = NULL;

	if (csv_file_noconv != NULL) {
		if ((csv_noconvPtr = fopen(csv_file_noconv,"w")) == NULL) {
			fprintf(stdout,"Could not write to %s\n",csv_file_noconv);
			return 1;
		}
		else if (options->verbose)
			g_fprintf(stdout,"Writing to CSV file %s\n",csv_file_noconv);
	}
	if (csv_file_conv != NULL) {
		if ((csv_convPtr = fopen(csv_file_conv,"w")) == NULL) {
			fprintf(stdout,"Could not write to %s\n",csv_file_conv);
			return 1;
		}
		else if (options->verbose)
			g_fprintf(stdout,"Writing to CSV file %s\n",csv_file_conv);
	}

	for (i =(zero_sum > 0.0 ? 0 : 1) ; i <= xi->general->n_interactions_trajectory ; i++) {

		//write it to outputfile... spe style
		if (spe_file_noconv != NULL) {
			filename = g_strdup_printf("%s_%i.spe",spe_file_noconv,i);
			if ((outPtr=fopen(filename,"w")) == NULL ) {
				fprintf(stdout,"Could not write to %s\n",filename);
				exit(1);
			}
			else if (options->verbose)
				g_fprintf(stdout,"Writing to SPE file %s\n",filename);

			fprintf(outPtr,"$SPEC_ID:\n\n");
				fprintf(outPtr,"$MCA_CAL:\n2\n");
				fprintf(outPtr,"%g %g\n\n", xi->detector->zero, xi->detector->gain);
			fprintf(outPtr,"$DATA:\n");
			fprintf(outPtr,"0\t%i\n",xi->detector->nchannels-1);
			for (j=0 ; j < xi->detector->nchannels ; j++) {
				fprintf(outPtr,"%g",ARRAY2D_FORTRAN(channels,i,j,xi->general->n_interactions_trajectory+1,xi->detector->nchannels));
				if ((j+1) % 8 == 0) {
					fprintf(outPtr,"\n");
				}
				else {
					fprintf(outPtr,"     ");
				}
			}
			fclose(outPtr);
			g_free(filename);
		}
		//convoluted spectrum
		if (spe_file_conv != NULL) {
			filename = g_strdup_printf("%s_%i.spe",spe_file_conv,i);
			if ((outPtr=fopen(filename,"w")) == NULL ) {
				fprintf(stdout,"Could not write to %s\n",filename);
				exit(1);
			}
			else if (options->verbose)
				g_fprintf(stdout,"Writing to SPE file %s\n",filename);

			fprintf(outPtr,"$SPEC_ID:\n\n");
			fprintf(outPtr,"$DATA:\n");
			fprintf(outPtr,"1\t%i\n",xi->detector->nchannels);
			for (j=0 ; j < xi->detector->nchannels ; j++) {
				fprintf(outPtr,"%g",channels_conv[i][j]);
				if ((j+1) % 8 == 0) {
					fprintf(outPtr,"\n");
				}
				else {
					fprintf(outPtr,"     ");
				}
			}
			fclose(outPtr);
			g_free(filename);
		}
	}

	xmi_free_input_F(&inputFPtr);

	//free stuff
	xmi_free_input(xi);
	free(xmim_input);
	free(spe_file);
	free(spe->mca_data);
	free(spe);
	free(weights_arr_quant);
	free(xmimsim_hdf5_escape_ratios);
	xmi_free_escape_ratios(escape_ratios_def);
	xmi_free_solid_angle(solid_angle_def);
	free(hdf5_file);
	xmi_deallocate(channels);
	xmi_deallocate(brute_history);
	xmi_deallocate(var_red_history);
	free(channels_conv);

	free(scale);
	free(k_exp);
	free(k_sim);
	free(l_exp);
	free(l_sim);
	free(energy_k);
	free(energy_l);
	free(k_sim_low);
	free(l_sim_low);
	free(k_sim_high);
	free(l_sim_high);
	free(weight_low);
	free(weight_high);

//	//Let's spawn xmim and perform calculation for the current input file
//	xmim_command = malloc(strlen("xmimsim --enable-opencl --enable-pile-up --spe-file xmim_opt ") + strlen(xmim_input) + 1);
//	if(xmim_command){
//		strcpy(xmim_command, "xmimsim --enable-opencl --enable-pile-up --spe-file xmim_opt ");
//		strcat(xmim_command, xmim_input);
//	}
//	printf("Command to run: %s\n",xmim_command);
//	fp = popen(xmim_command, "r");
//	while(fgets(str, 60, fp)!=NULL){
//		//xmim program still running
//	}
//	if(pclose(fp) == -1){
//		exit(0);
//	}



	return 0;

}
