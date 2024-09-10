#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <string>

// 计算并输出平均值到文件average.out
void averageToFile(const std::vector<std::string>& fileNames) {
    std::vector<std::ifstream> files;
    // 尝试打开所有文件
    for (const auto& name : fileNames) {
        files.emplace_back(name);
        if (!files.back()) {
            std::cerr << "无法打开文件 " << name << "\n";
            return;
        }
    }

    std::ofstream outFile("average.out");
    if (!outFile) {
        std::cerr << "无法创建输出文件average.out\n";
        return;
    }

    std::string line;
    while (true) {
        std::vector<double> values;
        double num;
        bool endOfFile = false;
        for (auto& file : files) {
            if (!std::getline(file, line)) {
                endOfFile = true;
                break;
            }
            std::istringstream iss(line);
            while (iss >> num) {
                values.push_back(num);
            }
        }
        if (endOfFile) break;

        int numFiles = fileNames.size();
        int numsPerFile = values.size() / numFiles;
        for (int i = 0; i < numsPerFile; ++i) {
            double sum = 0;
            for (int j = 0; j < numFiles; ++j) {
                sum += values[i + j * numsPerFile];
            }
            outFile << std::fixed << std::setprecision(15) << sum / numFiles << "\t";
        }
        outFile << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "用法：" << argv[0] << " file1 [file2 ...]\n";
        return 1;
    }
    std::vector<std::string> fileNames(&argv[1], &argv[argc]);
    averageToFile(fileNames);
    return 0;
}