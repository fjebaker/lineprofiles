DIST_DIR = .
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

LIB_EXT = a

all: _compile_for_xspec $(DIST_DIR)/kerr-transfer-functions.fits
	@echo "Successfully compiled"
	@echo "Use"
	@echo ""
	@echo "    lmod xsklineprofiles ./$(DIST_DIR)"
	@echo ""
	@echo "to load into XSPEC."

$(DIST_DIR)/kerr-transfer-functions.fits: kerr-transfer-functions.fits
	cp "$<" "$@"

kerr-transfer-functions.fits:
	curl -L \
		https://github.com/fjebaker/lineprofiles/releases/download/v0.1.0/kerr-transfer-functions-v0.1.0.zip \
		--output $(DIST_DIR)/kerr-transfer-functions.zip
	(cd $(DIST_DIR) && unzip kerr-transfer-functions.zip)

_compile_for_xspec: $(DIST_DIR)/libxsklineprofiles.$(LIB_EXT)
	(cd $(DIST_DIR) && \
		echo "initpackage xsklineprofiles lmodel.dat .\n exit" | xspec)
	rm $(DIST_DIR)/libxsklineprofiles.so
	$(SED_INPLACE) 's|-lXSFunctions|-lXSFunctions -l:libxsklineprofiles.$(LIB_EXT)|g' \
		$(DIST_DIR)/Makefile
	(cd $(DIST_DIR) && echo "hmake \n exit" | xspec)

.PHONY: clean
clean:
	rm -rf $(DIST_DIR)
	rm -rf zig-out
