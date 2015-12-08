SRCS = src/Utils.cpp src/WordCounter.cpp src/WordCount.cpp src/SequenceCounter.cpp src/SequenceCount.cpp src/GloveVLBL.cpp src/TrainAntonym.cpp src/TestAntonym.cpp
WC_OBJS = src/Utils.o src/WordCounter.o src/WordCount.o
SE_OBJS = src/Utils.o src/SequenceCounter.o src/SequenceCount.o
TR_OBJS = src/Utils.o src/TrainAntonym.o src/GloveVLBL.o 
TE_OBJS = src/Utils.o src/TestAntonym.o src/GloveVLBL.o 
OBJS = $(SRCS:.cpp=.o)
#CPPFLAGS= -fopenmp -pg -Og -DNDEBUG -std=c++11 -Wall 
CPPFLAGS= -DNDEBUG -fopenmp -fomit-frame-pointer -msse4.2 -mtune=native -march=native -std=c++11 -Wall -Ofast -funroll-loops -Wno-unused-result 
#CPPFLAGS= -fopenmp -fomit-frame-pointer -msse4.2 -mtune=native -march=native -std=c++11 -Wall -Ofast -funroll-loops -Wno-unused-result 
#CPPFLAGS= -fomit-frame-pointer -msse4.2 -mtune=native -march=native -std=c++11 -Wall -Ofast
#CPPFLAGS= -fopenmp -g -O2 -Wall -std=c++11
DEP = Makefile.depend
#LIBS = -L"/usr/local/lib" -lgzstream -lz -lboost_filesystem-mt -lboost_system-mt
LIBS = -lgzstream -lz -lm -lboost_filesystem -lboost_system
TARGETS = WordCount SequenceCount TrainAntonym TestAntonym


all: $(DEP) $(TARGETS)

WordCount: $(WC_OBJS)
	g++ $(CPPFLAGS) -o WordCount $(WC_OBJS) $(LIBS) 

SequenceCount: $(SE_OBJS)
	g++ $(CPPFLAGS) -o SequenceCount $(SE_OBJS) $(LIBS) 

TestAntonym: $(TE_OBJS)
	g++ $(CPPFLAGS) -o TestAntonym $(TE_OBJS) $(LIBS) 

TrainAntonym: $(TR_OBJS)
	g++ $(CPPFLAGS) -o TrainAntonym $(TR_OBJS) $(LIBS)


depend : $(DEP)
$(DEP): $(SRCS)
	g++ -MM -MG -std=c++11 $(SRCS)|sed -e 's/^\([^ ]\)/src\/\1/' > $(DEP)

clean:
	rm $(TARGETS) $(DEP) $(OBJS)

-include $(DEP)

