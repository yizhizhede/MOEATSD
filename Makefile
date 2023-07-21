#Makefile file

# compiler 
CC := g++ 

# option
FLAGS = -Wall -Iinclude # -no-multibyte-chars

# include director 
VPATH = include : src : bin

#
.PHONY : all clean 

# 
all : main moea 

#

# source file
src := $(filter %.c %.cxx %.cpp %.cu %.mpi,$(shell find .))

# dependent file
dependents := $(addsuffix .d,$(basename $(src)))
-include $(dependents)

# objective file 
objectives := $(addsuffix .o,$(basename $(src)))

#
main : $(objectives)
	$(CC) $(FLAGS) -Dpart_debug -g $^ -o bin/$@ -lm 
moea : 
	$(CC) $(FLAGS) -Dpart_debug -Dpart_release -O3 $(src) -o bin/$@ -lm
# 
%.o : %.c
	$(CC) -c $(FLAGS) -Dpart_debug $< -o $@
%.o : %.cpp
	$(CC) -c $(FLAGS) -Dpart_debug $< -o $@
%.o : %.cxx
	$(CC) -c $(FLAGS) -Dpart_debug $< -o $@
%.o : %.mpi
	$(CC) -c $(FLAGS) -Dpart_debug $< -o $@

# 
%.d : %.c
	@set -e; rm -f $@; \
	$(CC) -MM $(FLAGS) $< > $@.$$$$; \
	sed 's,\(.*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
%.d : %.cpp
	@set -e; rm -f $@; \
	$(CC) -MM $(FLAGS) $< > $@.$$$$; \
	sed 's,\(.*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
%.d : %.cxx
	@set -e; rm -f $@; \
	$(CC) -MM $(FLAGS) $< > $@.$$$$; \
	sed 's,\(.*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
%.d : %.mpi
	@set -e; rm -f $@; \
	$(CC) -MM $(FLAGS) $< > $@.$$$$; \
	sed 's,\(.*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
#
clean:
	-rm -f $(objectives) $(dependents) ./bin/*  
