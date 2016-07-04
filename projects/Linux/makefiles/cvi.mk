#INCLUDE
include make.inc


# COMPILER
PATH_OBJECTS 	= objects/$(COMPILER_TAG)/cvi
PATH_CPP 	= ../../src
PATH_HPP 	= ../../src

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/solver1d.o



# MAIN PROJECT
all: directories $(PATH_EXE)/OpenSMOKEpp_CVI$(STRING_NAME).sh

# CLEAN EVERYTHING
clean:
	-rm $(PATH_OBJECTS)/*.o $(PATH_EXE)/OpenSMOKEpp_CVI$(STRING_NAME).sh

# CREATE DIRECTORIES
directories: 	
	test -d $(PATH_EXE) || mkdir -p $(PATH_EXE)
	test -d $(PATH_OBJECTS) || mkdir -p $(PATH_OBJECTS)



# MAIN TARGET
$(PATH_EXE)/OpenSMOKEpp_CVI$(STRING_NAME).sh : $(OBJS)
	$(CXX) $(CXX_FLAGS) -o $(PATH_EXE)/OpenSMOKEpp_CVI$(STRING_NAME).sh $(OBJS) $(GLOBAL_LIB_INCLUDE) $(GLOBAL_LIBS)

# COMPILE FILES
$(PATH_OBJECTS)/solver1d.o : $(PATH_CPP)/solver1d.cpp
		 	         $(CCP_COMPILE) $(PATH_CPP)/solver1d.cpp -o $(PATH_OBJECTS)/solver1d.o


