/*
    レイトレーシングを利用して、OBJ形状をテキストでレンダリングするプログラム
    usage: ./ray-tracing <fname> <width> <height>
    memo: kd-treeを使えばさらに高速に処理できる
*/

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#define TIME_OUTPUT false

using namespace std;
using namespace Eigen;

bool ray_triangle_intersection(Vector3d const& a, Vector3d const& b,
                               Vector3d const& c, Vector3d const& p,
                               Vector3d const& d, Vector3d& u);
bool ray_mesh_intersection(vector<Vector3d> const& vertices,
                           vector<Vector3i> const& faces, Vector3d const& p,
                           Vector3d const& d, Vector3d& u, int& f_idx);
void ritter(vector<Vector3d> const& vertices, vector<Vector3i> const& faces,
            Vector3d& center, double& radius);
bool read(string const& fname, vector<Vector3d>& vertices,
          vector<Vector3i>& faces);
bool read_obj(string const& fname, vector<Vector3d>& vertices,
              vector<Vector3i>& faces);
bool read_off(string const& fname, vector<Vector3d>& vertices,
              vector<Vector3i>& faces);
void print_image(vector<vector<int>> const& image);

int main(int argc, char* argv[]) {
    // コマンドライン引数
    int width, height;
    string fname = "sample.obj";

    if (argc < 4) {
        cerr << "引数が足りません" << endl;
        cerr << "usage: " << argv[0] << " <fname> <width> <height>" << endl;
        return 1;
    }

    fname = argv[1];
    width = atoi(argv[2]);
    height = atoi(argv[3]);

    // メッシュ入力
    vector<Vector3d> vertices;
    vector<Vector3i> faces;

    if (!read(fname, vertices, faces)) {
        cerr << "read failed" << endl;
        return 1;
    }

    // 時間計測開始
    auto time0 = std::chrono::system_clock::now();

    // 最小包含円で正規化
    Vector3d center;
    double radius;
    ritter(vertices, faces, center, radius);
    for (Vector3d& v : vertices) v = (v - center) / radius;

    // 頂点法線計算
    vector<Vector3d> normals(vertices.size(), Vector3d::Zero());
    for (Vector3i f : faces) {
        Vector3d tmp = (vertices[f[1]] - vertices[f[0]])
                           .cross(vertices[f[2]] - vertices[f[0]]);
        normals[f[0]] += tmp;
        normals[f[1]] += tmp;
        normals[f[2]] += tmp;
    }
    for (Vector3d& n : normals) n.normalize();

    // レンダリング
    double ambient = 0.1, diffuse = 0.9;
    Vector3d light_vec = {1., 1., 1.};
    light_vec.normalize();
    vector<vector<int>> image(height,
                              vector<int>(width, 0xFFFFFF));  // 背景は白

    // レイトレーシング
    Vector3d p, d;
    Vector3d u;  // 重心座標
    int f_idx;   // 面番号
    unsigned char r, g, b;
    for (int h = 0; h < height; h++)
        for (int w = 0; w < width; w++) {
            // 色計算
            p = {-1. + (2. / width) * (w + 0.5),
                 -1. + (2. / height) * (height - h - 0.5), 2};
            d = {0, 0, -1};
            if (ray_mesh_intersection(vertices, faces, p, d, u, f_idx)) {
                // 交差あり
                Vector3d normal = Vector3d::Zero();
                normal += u[0] * normals[faces[f_idx][0]];
                normal += u[1] * normals[faces[f_idx][1]];
                normal += u[2] * normals[faces[f_idx][2]];
                normal.normalize();

                double l = ambient + diffuse * max(light_vec.dot(normal), 0.0);
                int c = clamp((int)(l * 256), 0, 255);
                r = g = b = c;
                image[h][w] = r << 16 | g << 8 | b;
            }
        }

    // 時間計測終了
    auto time1 = std::chrono::system_clock::now();
    using namespace std::chrono_literals;

    // 時間出力
    if (TIME_OUTPUT)
        cout << "time: " << (time1 - time0) / 1.0s << " sec" << endl;

    // 画像表示
    print_image(image);

    return 0;
}

// 最小包含円
void ritter(vector<Vector3d> const& vertices, vector<Vector3i> const& faces,
            Vector3d& center, double& radius) {
    // Ritters' method
    center = Vector3d::Zero();

    // 初期点を選択
    int p0 = 0;

    // p0から最も遠い点paを選ぶ
    int pa = p0;
    for (int i = 0; i < vertices.size(); i++)
        if ((vertices[i] - vertices[p0]).norm() >
            (vertices[pa] - vertices[p0]).norm())
            pa = i;

    // paから最も遠い点pbを選ぶ
    int pb = pa;
    for (int i = 0; i < vertices.size(); i++)
        if ((vertices[i] - vertices[pa]).norm() >
            (vertices[pb] - vertices[pa]).norm())
            pb = i;

    // 初期球を設定
    center = (vertices[pa] + vertices[pb]) / 2;
    radius = (vertices[pa] - center).norm();

    // 他の点について演算
    for (int i = 0; i < vertices.size(); i++) {
        // vertices[i] が S に含まれないときは球を更新
        if ((vertices[i] - center).norm() > radius) {
            double d = (vertices[i] - center).norm() - radius;
            center +=
                d / 2 * (vertices[i] - center) / (vertices[i] - center).norm();
            radius += d / 2;
        }
    }
}

// レイとメッシュの交点を求める
bool ray_mesh_intersection(vector<Vector3d> const& vertices,
                           vector<Vector3i> const& faces, Vector3d const& p,
                           Vector3d const& d, Vector3d& u, int& f_idx) {
    bool is_intersect = false;
    int min = 10;
    Vector3d tmp;
    for (int i = 0; i < faces.size(); i++) {
        if (ray_triangle_intersection(vertices[faces[i][0]],
                                      vertices[faces[i][1]],
                                      vertices[faces[i][2]], p, d, tmp)) {
            if ((vertices[faces[i][1]] - vertices[faces[i][0]])
                    .cross(vertices[faces[i][2]] - vertices[faces[i][0]])
                    .dot(d) < 0) {
                is_intersect = true;
                Vector3d intersection = vertices[faces[i][0]] * tmp[0] +
                                        vertices[faces[i][1]] * tmp[1] +
                                        vertices[faces[i][2]] * tmp[2];
                int l = (intersection - p).norm();
                if (l < min) {
                    min = l;
                    f_idx = i;
                    u = tmp;
                }
            }
        }
    }

    return is_intersect;
}

// レイと三角形の交点を求める
bool ray_triangle_intersection(Vector3d const& a, Vector3d const& b,
                               Vector3d const& c, Vector3d const& p,
                               Vector3d const& d, Vector3d& u) {
    // 重心座標
    Matrix3d mat;
    mat.col(0) = d;
    mat.col(1) = b - p;
    mat.col(2) = c - p;
    u[0] = mat.determinant();

    mat.col(1) = c - p;
    mat.col(2) = a - p;
    u[1] = mat.determinant();

    mat.col(1) = a - p;
    mat.col(2) = b - p;
    u[2] = mat.determinant();

    // 重心座標の和を1にする
    u /= u.sum();

    // 交差判定
    for (int i = 0; i < 3; i++)
        if (u[i] < 0) return false;
    return true;
}

// ファイル読み込み
bool read(string const& fname, vector<Vector3d>& vertices,
          vector<Vector3i>& faces) {
    string ext = filesystem::path(fname).extension().string();
    if (ext == ".obj")
        return read_obj(fname, vertices, faces);
    else if (ext == ".off")
        return read_off(fname, vertices, faces);
    return false;
}

// OBJファイルを読み込む関数
bool read_obj(string const& fname, vector<Vector3d>& vertices,
              vector<Vector3i>& faces) {
    ifstream ifs(fname);
    if (!ifs) {
        cerr << "file \"" << fname << "\" cannot be found!" << endl;
        return false;
    }

    // 頂点数と面数の取得
    int n = 0, f = 0;
    string line;
    while (getline(ifs, line)) {
        if (line.empty()) continue;
        switch (line[0]) {
            case 'v':
                n++;
                break;
            case 'f':
                f++;
                break;
        }
    };
    ifs.clear();

    vertices.resize(n);
    faces.resize(f);

    // 頂点情報と面情報の取得
    int cnt_n = 0, cnt_f = 0;
    string token;
    ifs.seekg(0);
    while (getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        istringstream iss(line);
        getline(iss, token, ' ');
        switch (line[0]) {
            case 'v':
                for (int i = 0; i < 3;) {
                    getline(iss, token, ' ');
                    if (token.empty()) continue;
                    vertices[cnt_n][i++] = stod(token);
                }
                cnt_n++;
                break;

            case 'f':
                for (int i = 0; i < 3;) {
                    getline(iss, token, ' ');
                    if (token.empty()) continue;
                    faces[cnt_f][i++] = stoi(token) - 1;
                }
                cnt_f++;
                break;

            default:
                continue;
                break;
        }
    }

    return true;
}

// OFFファイルを読み込む関数
bool read_off(string const& fname, vector<Vector3d>& vertices,
              vector<Vector3i>& faces) {
    ifstream ifs(fname);
    if (!ifs) {
        cerr << "file \"" << fname << "\" cannot be found!" << endl;
        return false;
    }

    // 頂点数と面数の取得
    int n = 0, f = 0;
    string line;
    string token;
    getline(ifs, line);
    getline(ifs, line);

    istringstream stream(line);
    getline(stream, token, ' ');
    n = stoi(token);
    getline(stream, token, ' ');
    f = stoi(token);

    vertices.resize(n);
    faces.resize(f);

    // 頂点情報の取得
    for (int i = 0; i < n; i++) {
        getline(ifs, line);
        istringstream iss(line);
        for (int j = 0; j < 3;) {
            getline(iss, token, ' ');
            if (token.empty()) continue;
            vertices[i][j++] = stod(token);
        }
    }

    for (int i = 0; i < f; i++) {
        getline(ifs, line);
        istringstream iss(line);
        getline(iss, token, ' ');
        for (int j = 0; j < 3;) {
            getline(iss, token, ' ');
            if (token.empty()) continue;
            faces[i][j++] = stoi(token);
        }
    }

    return true;
}

// 画像表示
void print_image(vector<vector<int>> const& image) {
    int height = image.size(), width = image[0].size(), px;

    for (int h = 0; h < height; h += 2) {
        for (int w = 0; w < width; w++) {
            // 上半分
            px = image[h][w];
            cout << "\e[38;2;" << ((px & 0xFF0000) >> 16) << ";"
                 << ((px & 0xFF00) >> 8) << ";" << (px & 0xFF) << "m";

            if (height % 2 == 0 || h != height - 1) {
                // 下半分
                px = image[h + 1][w];
                cout << "\e[48;2;" << ((px & 0xFF0000) >> 16) << ";"
                     << ((px & 0xFF00) >> 8) << ";" << (px & 0xFF) << "m";
            }

            // 描画
            cout << "\u2580";
        }
        cout << "\e[0m" << endl;
    }
}
