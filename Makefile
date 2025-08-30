.PHONY: build clean

build:
	mkdir -p build
	cd build && cmake .. -GNinja -DCXX=clang++ -DCC=clang
	cd build && cmake --build .

clean:
	rm -rf build
