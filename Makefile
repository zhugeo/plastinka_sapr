.PHONY: build clean

configure:
	mkdir -p build
	cd build && cmake .. -GNinja -DCXX=clang++ -DCC=clang

build: configure
	cd build && cmake --build .

clean:
	rm -rf build
