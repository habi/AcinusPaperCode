all: dtSkel

COMMON_DIR=../Common
UTILS_DIR=../Utils/include

common.o: ${COMMON_DIR}/common.cpp ${COMMON_DIR}/common.h 
	g++ -c ${COMMON_DIR}/common.cpp -o common.o -I${COMMON_DIR}/

dtSkel: dtSkel.cpp getDT.h common.o
	g++ dtSkel.cpp common.o -o dtSkel -I${COMMON_DIR}/

thin: thinning-main.cpp thinning.h common.o ${UTILS_DIR}/openFile.h
	g++ thinning-main.cpp common.o -o thin -I${UTILS_DIR} -I${COMMON_DIR}

clean: 
	rm -f dtSkel thin *~ *.o
