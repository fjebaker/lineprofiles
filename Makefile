DIST_DIR = kline-xspec

all:
	rm -rf $(DIST_DIR)
	zig build --release=fast -Dtarget=x86_64-linux-gnu xspec
	cp -r ./xspec $(DIST_DIR)
	cp ./zig-out/lib/libxsklineprofiles.a $(DIST_DIR)
	cp ./kerr-transfer-functions.fits $(DIST_DIR)
	make -C $(DIST_DIR) -f MakeForXSPEC

.PHONY: clean
clean:
	rm -rf $(DIST_DIR)
