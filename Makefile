run:
	cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
	cmake --build build-debug
	./build-debug/transmission_line_waves

run-release:
	cmake -S . -B build
	cmake --build build 
	./build/transmission_line_waves
