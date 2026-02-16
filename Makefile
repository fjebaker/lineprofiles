DIST_DIR = kline-xspec
TARGET = $(shell uname)
ARCH = $(shell uname -m)

ifeq ($(TARGET),Linux)
	OSTARGET=$(ARCH)-linux-gnu
	SHARED_EXT := so
	SED_INPLACE = sed -i
else
ifeq ($(TARGET),Darwin)
	OSTARGET=$(ARCH)-macos-none
	SHARED_EXT := dylib
	SED_INPLACE = sed -i ''
endif
endif

LIB_PATH := $(abspath $(DIST_DIR))

# TODO: why is it this strange file extension?
DOWNLOAD_TARGET = $(OSTARGET)-lineprofiles.tar.gz.zip

all: _compile_for_xspec data

data: $(DIST_DIR)/kerr-transfer-functions.fits

$(DIST_DIR)/kerr-transfer-functions.fits:
	curl -L \
		https://github.com/fjebaker/lineprofiles/releases/download/v0.1.0/kerr-transfer-functions-v0.1.0.zip \
		--output $@

_compile_for_xspec: $(DIST_DIR)/libxsklineprofiles.a
	# delete everything that is not needed
	rm -f \
		 $(DIST_DIR)/lpack_klineprofiles.cxx \
		 $(DIST_DIR)/lpack_klineprofiles.o \
		 $(DIST_DIR)/Makefile \
		 $(DIST_DIR)/pkgIndex.tcl \
		 $(DIST_DIR)/x86_64-linux-gnu-lineprofiles.tar.gz.zip \
		 $(DIST_DIR)/klineprofilesFunctionMap.cxx \
		 $(DIST_DIR)/klineprofilesFunctionMap.h \
		 $(DIST_DIR)/klineprofilesFunctionMap.o
	(cd $(DIST_DIR) && \
		echo "initpackage xsklineprofiles lmodel.dat .\n exit" | xspec)
	rm $(DIST_DIR)/libxsklineprofiles.$(SHARED_EXT)
	$(SED_INPLACE) 's|-lXSFunctions|-lXSFunctions -L$(LIB_PATH) -Wl,-rpath,"$(LIB_PATH)" -l:libxsklineprofiles.a|g' \
		$(DIST_DIR)/Makefile
	(cd $(DIST_DIR) && echo "hmake \n exit" | xspec)

$(DIST_DIR)/libxsklineprofiles.a: $(DOWNLOAD_TARGET)
	rm -rf $(DIST_DIR)
	mkdir -p $(DIST_DIR)
	cp $(DOWNLOAD_TARGET) $(DIST_DIR)
	(cd $(DIST_DIR) && unzip $(DOWNLOAD_TARGET) && rm -f ./$(DOWNLOAD_TARGET))

$(DOWNLOAD_TARGET):
	curl -L \
		https://github.com/fjebaker/lineprofiles/releases/download/v0.1.2/$(DOWNLOAD_TARGET) \
		--output $(DOWNLOAD_TARGET)

.PHONY: xspec-test
xspec-test:
	rm -rf build
	cp -r xspec build
	cp zig-out/lib/libxsklineprofiles.a build
	(cd build && echo "hmake \n exit \n" | xspec)
	cp ./kerr-transfer-functions.fits build

.PHONY: clean
clean:
	rm -rf $(DIST_DIR)
