CXX      := clang++
CXXFLAGS := -std=c++17 -O3 -Icore/eigen
TARGET   := quadra_example
SRC      := main.cpp
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.o a.out

.PHONY: all run clean
