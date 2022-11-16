out/main: *.cpp *.h
	mkdir -p out
	g++ -std=c++11 -Wall -Wextra -Wfatal-errors -O3 \
 		-Wpedantic -pedantic-errors \
		main.cpp \
		-o out/main
	
clean:
	rm -rf out
