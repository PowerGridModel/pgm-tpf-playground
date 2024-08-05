#include <common.hpp>

int main() {
    std::vector<std::string> names = {"PGM", "TPF", "Hackathon", "2024"};
    for (const auto& name : names) {
        std::cout << name << std::endl;
    }
    return 0;
}