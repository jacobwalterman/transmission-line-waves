run:
	cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
	cmake --build build-debug
	./build-debug/coupled_piezo_heat

run-release:
	cmake -S . -B build
	cmake --build build 
	./build/coupled_piezo_heat
