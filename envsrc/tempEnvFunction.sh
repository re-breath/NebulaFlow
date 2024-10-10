# 该库使用来放置临时函数，多为针对性极强的函数，泛用性一般没有


del_caf2xyz_F(){
# 删除caf2文件中的F原子,用于处理xyz文件
filename=${1:-"dump.xyz"}
file2name=${2:-"del_${filename}"}
> caf2_tmpfile.cpp
cat >> caf2_tmpfile.cpp <<EOF
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>

using namespace std;

pair<string, string> get_needReplace(const string& filename) {
    ifstream file(filename);
    string firstLine;
    getline(file, firstLine);  // 读取第一行
    int needReplace = stoi(firstLine);
    string replace2 = to_string(needReplace / 3);
    string needReplaceStr = to_string(needReplace);
    file.close();
    return make_pair(needReplaceStr, replace2);
}

void del_F(const string& filename) {
    ifstream input_file(filename);
    ofstream output_file("$file2name");

    pair<string, string> result = get_needReplace(filename);
    string needReplace = result.first;
    string replace2 = result.second;

    string line;
    while (getline(input_file, line)) {
        if (line.empty()) {
            continue;  // 跳过空行
        }

        stringstream ss(line);
        string token;
        vector<string> tokens;
        while (getline(ss, token, ' ')) {
            tokens.push_back(token);
        }

        if (tokens[0] == needReplace) {
            tokens[0] = replace2;
        } else if (tokens[0] == "F") {
            continue;
        }

        for (size_t i = 0; i < tokens.size(); ++i) {
            output_file << tokens[i];
            if (i != tokens.size() - 1) {
                output_file << " ";
            }
        }
        output_file << "\n";
    }

    input_file.close();
    output_file.close();
}

int main() {
    string filename = "$filename";

    auto start_time = std::chrono::high_resolution_clock::now();
    del_F(filename);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

    cout << "Using Time: " << elapsed_time << " ms" << endl;

    return 0;
}
EOF

g++  -O3 caf2_tmpfile.cpp -o caf2_tmpfile
./caf2_tmpfile
rm caf2_tmpfile.cpp caf2_tmpfile
}


verify_deform(){
#临时函数，用来验算拉伸曲线的结果
    verify_gpumd_result
    cd verify_*
    sed -i  "/dump_thermo/c\dump_thermo  100" run.in
    sed -i  "/run/c\run 1000000" run.in
    gpumdstart
}

