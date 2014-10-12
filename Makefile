# gtrieScanner: quick discovery of network motifs
#
# Pedro Ribeiro - CRACS & INESC-TEC, DCC/FCUP
# ------------------------------------------------
# Makefile
# ------------------------------------------------

# Executable name
EXEC_NAME=gtrieScanner

# Compiler to use
CC=g++

# Default Options to use
# CFLAGS= -Wall -Wno-write-strings -O0 -g
CFLAGS= -Wno-write-strings -O3
#CFLAGSFINAL= -O3
CLIBS=

#NAUTYSRC= nauty.c nautil.c naugraph.c

# Source files
SRC =                   \
	CmdLine.cpp	\
	GraphMatrix.cpp	\
	GraphUtils.cpp	\
	Esu.cpp		\
	Conditions.cpp	\
	Isomorphism.cpp	\
	Error.cpp	\
	GTrie.cpp	\
	Timer.cpp	\
	Random.cpp	\
	GraphTree.cpp	\
	nauty/nauty.c	\
	nauty/nautil.c	\
	nauty/naugraph.c	\
	main.cpp

OBJ =  ${SRC:.cpp=.o}

#------------------------------------------------------------

all: ${EXEC_NAME}

${EXEC_NAME}: ${OBJ}
	${CC} ${CFLAGS} ${CLIBS} -o ${EXEC_NAME} ${OBJ}

%.o: %.cpp
	${CC} ${CFLAGS} -c -o $@ $+

clean:
	rm ${EXEC_NAME} *.o *~ *# -rf
