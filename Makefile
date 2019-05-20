CFLAGS = -fdiagnostics-color=always -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include -I/usr/local/include/xmimsim -I/usr/local/include/xraylib -I/usr/include/libxml2 -I/usr/local/include -g -O0 -Wall
CC=gcc
LIBS = -lm -lglib-2.0 -lxrl -lxmimsim -lglib-2.0 -lm
OBJ = xmim_auto.o 

xmim-auto: $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@  $^ $(LIBS)
