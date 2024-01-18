#include <iostream>
#include <vector>

int main() {
    std::vector<int> myVector = {1, 2, 3, 4, 5};

    // 删除第一个元素
    if (!myVector.empty()) {
        myVector.erase(myVector.begin());
    }

    // 打印更新后的std::vector
    for (const auto &element : myVector) {
        std::cout << element << " ";
    }

    return 0;
}