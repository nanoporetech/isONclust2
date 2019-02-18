buildDir ?= build
releaseType ?= Release

.PHONY: all isONclust2
all: isONclust2
isONclust2: ${buildDir}/bin/isONclust2

${buildDir}:
	mkdir ${buildDir}

.PHONY: clean
clean:
	rm -rf ${buildDir}

.PHONY: docs
docs:
	doxygen Doxyfile

${buildDir}/bin/isONclust2: ${buildDir}
	cd ${buildDir} && \
	cmake .. -DCMAKE_BUILD_TYPE=${releaseType} && \
	make -j

.PHONY: test
test: ${buildDir}/bin/isONclust2
	@${buildDir}/bin/test_isONclust2
