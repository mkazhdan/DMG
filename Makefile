CLIENT_TARGET=Client
SERVER_TARGET=Server
BIG_IMAGE_PROCESS_TARGET=BigImageProcess
SOURCE= \
	Util/BaseMultiStreamIO.cpp \
	Util/CmdLineParser.cpp \
	Util/GridStream.cpp \
	Util/Half/half.cpp \
	Util/MultiStreamIO.cpp \
	Util/Socket.cpp \
	Util/Time.cpp \
	Util/XPlatform.cpp

CLIENT_SOURCE=$(SOURCE) ClientSocket/ClientSocket.cpp
SERVER_SOURCE=$(SOURCE) ServerSocket/ServerSocket.cpp
BIG_IMAGE_PROCESS_SOURCE=$(SOURCE) BigImageProcess/BigImageProcess.cpp

COMPILER ?= gcc
#COMPILER ?= clang

ifeq ($(COMPILER),gcc)
	CFLAGS += -fpermissive -fopenmp -Wno-deprecated -msse2 --std=c++11 -Wno-unused-result -Wno-dangling-else -Wno-pointer-bool-conversion
	LFLAGS += -lgomp -lz -lpng -ltiff -ljpeg -lboost_thread -lboost_system -lpthread
else
	CFLAGS += -fpermissive -Wno-deprecated -msse2 --std=c++11 -Wno-unused-result -Wno-dangling-else -Wno-pointer-bool-conversion
	LFLAGS += -lz -lpng -ltiff -ljpeg -lboost_thread -lboost_system -lpthread
endif

CFLAGS_DEBUG = -DDEBUG -g3
LFLAGS_DEBUG =

CFLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -ffast-math
LFLAGS_RELEASE = -O3 

SRC = ./
BIN = Bin/
BIN_O = ./
INCLUDE = /usr/include/ -I.

ifeq ($(COMPILER),gcc)
	CC=gcc
	CXX=g++
else
	CC=clang
	CXX=clang++
endif

MD=mkdir


##CLIENT_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(CLIENT_SOURCE))))
##SERVER_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(SERVER_SOURCE))))
CLIENT_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(CLIENT_SOURCE))))
SERVER_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(SERVER_SOURCE))))
BIG_IMAGE_PROCESS_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(BIG_IMAGE_PROCESS_SOURCE))))

all: CFLAGS += $(CFLAGS_RELEASE)
all: LFLAGS += $(LFLAGS_RELEASE)
all: $(BIN)$(CLIENT_TARGET)
all: $(BIN)$(SERVER_TARGET)
all: $(BIN)$(BIG_IMAGE_PROCESS_TARGET)

debug: CFLAGS += $(CFLAGS_DEBUG)
debug: LFLAGS += $(LFLAGS_DEBUG)
debug: $(BIN)$(CLIENT_TARGET)
debug: $(BIN)$(SERVER_TARGET)
debug: $(BIN)$(BIG_IMAGE_PROCESS_TARGET)

clean:
	rm -f $(BIN)$(CLIENT_TARGET)
	rm -f $(BIN)$(SERVER_TARGET)
	rm -f $(BIN)$(BIG_IMAGE_PROCESS_TARGET)
	rm -f $(CLIENT_OBJECTS)
	rm -f $(SERVER_OBJECTS)
	rm -f $(BIG_IMAGE_PROCESS_OBJECTS)

$(BIN):
	mkdir -p $(BIN)

$(BIN)$(CLIENT_TARGET): $(CLIENT_OBJECTS)
	mkdir -p $(BIN)
	$(CXX) -o $@ $(CLIENT_OBJECTS) $(LFLAGS)

$(BIN)$(SERVER_TARGET): $(SERVER_OBJECTS)
	mkdir -p $(BIN)
	$(CXX) -o $@ $(SERVER_OBJECTS) $(LFLAGS)

$(BIN)$(BIG_IMAGE_PROCESS_TARGET): $(BIG_IMAGE_PROCESS_OBJECTS)
	mkdir -p $(BIN)
	$(CXX) -o $@ $(BIG_IMAGE_PROCESS_OBJECTS) $(LFLAGS)

$(BIN_O)%.o: $(SRC)%.c
	$(CC) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

$(BIN_O)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

