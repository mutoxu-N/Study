
/*
    OBJファイルを読み込んで、メッシュ分割を行うプログラム
    usage: ./split-mesh <input_fname> <output_fname> <spd> <median_cut>
*/

#include <Spectra/MatOP/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#define TIME_OUTPUT true

using namespace std;
using Vector3d = Eigen::Vector3d;
using VectorXd = Eigen::VectorXd;
using Vector3i = Eigen::Vector3i;
using Triplet = Eigen::Triplet<double>;
using SparseMatrix = Eigen::SparseMatrix<double>;
using SparseVector = Eigen::SparseVector<double>;
using SparseSymMatProd = Spectra::SparseSymMatProd<double>;
using SymEigsSolver = Spectra::SymEigsSolver<SparseSymMatProd>;

bool read_obj(string const& fname, vector<Vector3d>& vertices,
              vector<Vector3i>& faces);
void write_obj(string const& fname, vector<Vector3d> const& vertices,
               vector<Vector3i> const& faces, vector<Vector3d> const& normals,
               vector<char> const& group);

// Halfedge
int prev(int k) { return k % 3 == 0 ? k + 2 : k - 1; }

int main(int argc, char* argv[]) {
    // コマンドライン引数
    int width, height;
    int spd = 20;
    bool median_cut = false;
    string fname = "sample.obj";
    string fname_out = "sample_out.obj";

    if (argc > 1) {
        fname = argv[1];
    }
    if (argc > 2) {
        fname_out = argv[2];
    }
    if (argc > 3) {
        spd = stoi(argv[3]);
    }
    if (argc > 4) {
        if (stoi(argv[4]) != 0) median_cut = true;
    }

    // メッシュ入力
    vector<Vector3d> vertices;
    vector<Vector3i> faces;

    if (!read_obj(fname, vertices, faces)) {
        cerr << "read failed" << endl;
        return 1;
    }

    // 時間計測開始
    auto time0 = std::chrono::system_clock::now();

    // 面の法線ベクトル
    vector<Vector3d> face_normals(faces.size());
    for (int i = 0; i < faces.size(); i++) {
        Vector3i face = faces[i];
        Vector3d v0 = vertices[face[0]];
        Vector3d v1 = vertices[face[1]];
        Vector3d v2 = vertices[face[2]];
        Vector3d normal = (v1 - v0).cross(v2 - v0);
        face_normals[i] = normal.normalized();
    }

    // ハーフエッジデータ構造の構築
    vector<int> opps = vector(3 * faces.size(), -1);  // 半辺数
    vector<int> outs = vector(vertices.size(), -1);   // outgoing edge

    // outgoing edgeリストの構築
    int t;
    for (int i = 0; i < faces.size(); i++) {
        // 各面の各頂点に対して, h_outを書き込む
        t = 3 * i;
        outs[faces[i][0]] = t + 2;
        outs[faces[i][1]] = t;
        outs[faces[i][2]] = t + 1;
    }

    // opposite edgeリストの構築
    map<pair<int, int>, int> he_map;
    int src, tgt, k;
    t = 3 * faces.size();
    for (int i = 0; i < t; i++) {
        src = faces[i / 3][(i + 1) % 3];
        tgt = faces[i / 3][(i + 2) % 3];

        if (he_map.count(pair(tgt, src)) > 0) {
            // マップ内に逆向き半辺が存在する
            k = he_map[pair(tgt, src)];
            opps[k] = i;
            opps[i] = k;
            he_map.erase(pair(tgt, src));

        } else {
            // マップ内に逆向き半辺が無い
            he_map[pair(src, tgt)] = i;
        }
    }

    // 各境界辺の始点に対する h_out をその境界辺に設定する
    for (auto m = he_map.begin(); m != he_map.end(); m++)
        outs[m->first.first] = m->second;

    // グラフラプラシアンの構築
    // 同時に頂点の法線ベクトルも計算する
    vector<Vector3d> normals(vertices.size());
    SparseMatrix laplacian(vertices.size(), vertices.size());
    vector<Triplet> triplets;
    int ed, he, cnt;
    vector<int> diagonal(vertices.size(), 0);
    for (int he = 0; he < 3 * faces.size(); he++) {
        src = faces[he / 3][(he + 1) % 3];
        tgt = faces[he / 3][(he + 2) % 3];
        triplets.push_back(Triplet(src, tgt, -1));
        diagonal[src]++;
    }
    for (int i = 0; i < vertices.size(); i++) {
        triplets.push_back(Triplet(i, i, diagonal[i]));
    }
    laplacian.setFromTriplets(triplets.begin(), triplets.end());

    // Fiedlerベクトル
    SparseSymMatProd op(laplacian);
    SymEigsSolver solver(op, 2, spd);
    VectorXd fiedler;
    solver.init();
    solver.compute(Spectra::SortRule::SmallestAlge);

    if (solver.info() == Spectra::CompInfo::Successful) {
        // 2番目に小さな固有値の固有ベクトル
        cout << "2番目に小さな固有値: " << solver.eigenvalues()[0] << endl;
        cout << "最も小さな固有値: " << solver.eigenvalues()[1] << endl;
        fiedler = solver.eigenvectors().col(0);

    } else {
        cout << "固有値計算に失敗しました: ";
        switch (solver.info()) {
            case Spectra::CompInfo::NotComputed:
                cout << "Not Computed" << endl;
                break;
            case Spectra::CompInfo::NotConverging:
                cout << "NotConverging" << endl;
                break;
            case Spectra::CompInfo::NumericalIssue:
                cout << "NumericalIssue" << endl;
                break;
        }
        return 1;
    }

    // 分割
    vector<char> group(vertices.size());
    double threshold;
    if (median_cut) {
        vector<double> t(vertices.size());
        for (int i = 0; i < vertices.size(); i++) t[i] = fiedler[i];
        sort(t.begin(), t.end());
        if (vertices.size() % 2 == 0) {
            threshold =
                (t[vertices.size() / 2] + t[vertices.size() / 2 - 1]) / 2;

        } else {
            threshold = t[vertices.size() / 2];
        }

    } else {
        threshold = 0;
    }

    cout << "threshold: " << threshold << endl;

    int red = 0, blue = 0;
    for (int i = 0; i < vertices.size(); i++) {
        if (fiedler[i] > threshold) {
            group[i] = 1;
            red++;
        } else {
            group[i] = -1;
            blue++;
        }
    }

    // 時間計測終了
    auto time1 = std::chrono::system_clock::now();

    // カットの重み
    int cut_weight = 0;
    for (int i = 0; i < 3 * faces.size(); i++) {
        src = faces[i / 3][(i + 1) % 3];
        tgt = faces[i / 3][(i + 2) % 3];
        if (group[src] < group[tgt]) {  // 青→赤の反辺数をカウント
            cut_weight++;
        }
    }

    // 連結性の確認
    stack<int> s;
    vector<bool> done(vertices.size(), false);
    int color = group[0];
    int next;

    // 頂点0の色で深さ優先探索
    s.push(0);
    done[0] = true;
    while (s.size() > 0) {
        int v = s.top();
        s.pop();
        ed = outs[v];
        he = ed;

        // 頂点vの半辺を1週する
        do {
            tgt = faces[he / 3][(he + 2) % 3];
            if (!done[tgt]) {
                if (group[tgt] == color) {
                    s.push(tgt);
                    done[tgt] = true;

                } else {
                    // 異なる色の頂点を保存しておく
                    next = tgt;
                }
            }
            he = opps[prev(he)];
        } while (he >= 0 && he != ed);
    }

    // 頂点0と異なる色で深さ優先探索
    s.push(next);
    done[next] = true;
    while (s.size() > 0) {
        int v = s.top();
        s.pop();
        ed = outs[v];
        he = ed;

        // 頂点vの半辺を1週する
        do {
            tgt = faces[he / 3][(he + 2) % 3];
            if (!done[tgt] && group[tgt] != color) {
                s.push(tgt);
                done[tgt] = true;

            } else {
                next = tgt;
            }
            he = opps[prev(he)];
        } while (he >= 0 && he != ed);
    }

    // 各領域が連結ならdone内にfalseが存在しない
    cnt = 0;
    for (int i = 0; i < vertices.size(); i++) {
        if (!done[i]) cnt++;
    }

    // 出力
    cout << "vertices: " << vertices.size() << endl;
    cout << "faces: " << faces.size() << endl;
    cout << "red: " << red << endl;
    cout << "blue: " << blue << endl;
    cout << "cut_weight: " << cut_weight << endl;
    cout << "connectivity: " << (cnt == 0 ? "true" : "false") << endl;

    write_obj(fname_out, vertices, faces, normals, group);

    // 時間出力
    using namespace std::chrono_literals;
    if (TIME_OUTPUT)
        cout << "time: " << (time1 - time0) / 1.0s << " sec" << endl;
    return 0;
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

void write_obj(string const& fname, vector<Vector3d> const& vertices,
               vector<Vector3i> const& faces, vector<Vector3d> const& normals,
               vector<char> const& group) {
    ofstream ofs(fname);
    if (!ofs) {
        cerr << "file \"" << fname << "\" cannot be created!" << endl;
        return;
    }

    // 頂点書き込み
    for (int i = 0; i < vertices.size(); i++) {
        ofs << "v " << vertices[i].transpose();
        if (group[i] > 0) {
            ofs << " " << 1.0 << " " << 0.0 << " " << 0.0 << endl;
        } else {
            ofs << " " << 0.0 << " " << 0.0 << " " << 1.0 << endl;
        }
        ofs << "vn " << normals[i].transpose() << endl;
    }

    // 面書き込み
    for (int i = 0; i < faces.size(); i++) {
        ofs << "f " << faces[i][0] + 1 << " " << faces[i][1] + 1 << " "
            << faces[i][2] + 1 << endl;
    }

    ofs.close();
}
