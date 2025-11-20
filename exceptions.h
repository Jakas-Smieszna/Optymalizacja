#include <stdexcept>
#include <string>

class ToManyCalls : public std::exception {
private:
    string info;
public:
    ToManyCalls() {
        info = "Przekroczono limit wywolan funkcji celu.\n";
    }
    ToManyCalls(string info0) {
        info = ("Przekroczono limit wywolan funkcji celu:\n" + info0 + "\n");
    }
    const char* what() const noexcept override {
        return info.c_str();
    }
};

class BadArguments : public std::exception {
private:
    string info;
public:
    BadArguments() {
        info = "Otrzymane argumenty nie spelniaja wymogow algorytmu.\n";
    }
    BadArguments(string info0) {
        info = ("Otrzymane argumenty nie spelniaja wymogow algorytmu:\n" + info0 + "\n");
    }
    const char* what() const noexcept override {
        return info.c_str();
    }
};