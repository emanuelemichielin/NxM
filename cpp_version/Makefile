BLAS_INC?=.
BLAS_LIB?=.

CXX= g++-7
CXXFLAGS=-O3 -g -march=native -std=c++14 -std=gnu++14 -I${BLAS_INC} 
LDFLAGS=-g -L${BLAS_LIB} 
LDLIBS=-lopenblas -lfftw3 -llapacke -llapack `root-config --cflags --libs`
LD=${CXX}
ROOTFLAGS = -I${ROOTSYS}/include



output: main.o utils.o NoisePSDGenerator.o TemplateGeneratorNxM.o OptimalFilterNxM.o OFManagerNxM.o Driver.o
	${LD} ${LDFLAGS} -o reco main.o utils.o NoisePSDGenerator.o TemplateGeneratorNxM.o OptimalFilterNxM.o OFManagerNxM.o Driver.o ${LDLIBS}

main.o: main.cc header.h
	${CXX} ${CXXFLAGS} ${ROOTFLAGS} -c -o main.o main.cc

utils.o: utils.cc header.h
	${CXX} ${CXXFLAGS} ${ROOTFLAGS} -c -o utils.o utils.cc

NoisePSDGenerator.o: NoisePSDGenerator.cc header.h NoisePSDGenerator.h
	${CXX} ${CXXFLAGS} ${ROOTFLAGS} -c -o NoisePSDGenerator.o NoisePSDGenerator.cc

TemplateGeneratorNxM.o: TemplateGeneratorNxM.cc header.h TemplateGeneratorNxM.h
	${CXX} ${CXXFLAGS} ${ROOTFLAGS} -c -o TemplateGeneratorNxM.o TemplateGeneratorNxM.cc

OptimalFilterNxM.o: OptimalFilterNxM.cc header.h OptimalFilterNxM.h
	${CXX} ${CXXFLAGS} ${ROOTFLAGS} -c -o OptimalFilterNxM.o OptimalFilterNxM.cc

OFManagerNxM.o: OFManagerNxM.cc header.h OFManagerNxM.h
	${CXX} ${CXXFLAGS} ${ROOTFLAGS} -c -o OFManagerNxM.o OFManagerNxM.cc

Driver.o: Driver.cc header.h
	${CXX} ${CXXFLAGS} ${ROOTFLAGS} -c -o Driver.o Driver.cc

clean:
	\rm -f reco main.o NoisePSDGenerator.o TemplateGeneratorNxM.o OptimalFilterNxM.o OFManagerNxM.o utils.o Driver.o

clen:
	rm core.*
