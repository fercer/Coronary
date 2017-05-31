CCOMP=g++
CLINK=g++
FLAGS=-lm -DNDEBUG
OPENCV_FLAGS=`pkg-config --libs --cflags opencv`
FFTW3L_FLAGS=`pkg-config --libs --cflags fftw3l`
SDIR=./src
ODIR=./obj
PROG=./coronary

#Fuentes:
SRCS:=$(wildcard $(SDIR)/*.cpp)
HDRS:=$(wildcard $(SDIR)/*.h)
BINARIES:=$(PROG).bin

#Objetos:
OBJS:=$(patsubst $(SDIR)/%, $(ODIR)/%, $(patsubst %.cpp, %.o, $(SRCS)))
COMP_OBJS:=$(patsubst $(ODIR)/%, ./%, $(OBJS))

.PHONY: make
make: $(OBJS)
	$(CLINK) $(OBJS) $(HDRS) -o $(PROG).bin $(OPENCV_FLAGS) $(FFTW3L_FLAGS) $(FLAGS)

$(OBJS): $(SRCS) $(ODIR)/
	$(CCOMP) $(SRCS) -c $(OPENCV_FLAGS) $(FFTW3L_FLAGS) $(FLAGS)
	mv --target-directory=$(ODIR)/ $(COMP_OBJS)

$(ODIR):
	mkdir $(ODIR)

.PHONY: clean
clean:
	rm $(BINARIES) -f
	rm $(ODIR) -R -f

.PHONY: debug
debug:
	@echo "\n###################### DEBUG MODE ON ########################################\n"	
	make -f makefile FLAGS="-g -lm"
