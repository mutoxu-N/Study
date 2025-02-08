/*
    アンカー付き長方形の面積を最大化する遺伝的アルゴリズム

    usage: ./genetics <fname> <n>
    DATA_DIR の fname を読み込む 頂点数は $n$
*/

#include <bits/stdc++.h>
#include <direct.h>

#define E 10e-8
#define DIV 128
#define CHILDS 100
#define GENERATIONS 1000

using namespace std;
using P = array<double, 2>;
using V = array<double, 2>;
using R = array<P, 2>;

bool is_overlap(R const& a, R const& b);
bool is_overlap(R const& r, P const& p);
double rect_area(R const& r) {
    double dx = r[1][0] - r[0][0];
    if (abs(dx) < E) dx = 0;
    double dy = r[1][1] - r[0][1];
    if (abs(dy) < E) dy = 0;
    return dx * dy;
}
double area_check(vector<P> const& ps, vector<R> const& ans);
bool read(string const& fname, vector<P>& ps, vector<R>& rs);
array<double, 2> crossover(vector<P> const& ps, vector<R>& rs, mt19937& eng,
                           uniform_real_distribution<double>& dist,
                           vector<V>& gene1, vector<V>& gene2);
double simulate(vector<P> const& ps, vector<R>& rs, vector<V> const& gene);

int main(int argc, char* argv[]) {
    string const DATA_DIR = "../data/";
    string fname = "sample.txt";
    int n = 100;

    if (argc > 1) fname = argv[1];
    if (argc > 2) n = stoi(argv[2]);

    string const OUT_DIR = "./results/";
    string const OUT_FNAME = OUT_DIR + fname;
    _mkdir(OUT_DIR.c_str());

    vector<P> ps(n);
    vector<R> rs(n);
    assert(read(fname, ps, rs));

    // 時間計測開始
    auto time0 = std::chrono::system_clock::now();

    // アルゴリズム
    vector<V> gene1(n), gene2(n);

    mt19937 eng(0);
    uniform_real_distribution<double> dist(0, 1);

    for (int i = 0; i < n; i++) {
        gene1[i] = {dist(eng), dist(eng)};
        gene2[i] = {dist(eng), dist(eng)};
    }

    double area = 0;
    vector<V> best;
    for (int gen = 0; GENERATIONS - (gen++);) {
        array<double, 2> results = crossover(ps, rs, eng, dist, gene1, gene2);
        double t = results[0];
        if (area < t) {
            area = t;
            best = gene1;
        }
    }

    // 時間計測終了
    auto time1 = std::chrono::system_clock::now();

    double calc_area = simulate(ps, rs, best);
    double check_area = area_check(ps, rs);
    assert(abs(calc_area - check_area) < E);
    assert(abs(area - check_area) < E);

    cout << area << endl;

    // 実行時間
    using namespace std::chrono_literals;
    cout << (time1 - time0) / 1.0ms << "ms" << endl;

    // 出力
    cout << OUT_FNAME << endl;
    ofstream ofs(OUT_FNAME);
    for (R r : rs) {
        ofs << r[0][0] << " " << r[0][1] << " " << r[1][0] << " " << r[1][1]
            << endl;
    }
    ofs.close();

    return 0;
}

bool read(string const& fname, vector<P>& ps, vector<R>& rs) {
    ifstream ifs(fname);
    if (!ifs) {
        cerr << "file \"" << fname << "\" cannot be found!" << endl;
        return false;
    }

    int n = ps.size();
    assert(rs.size() == n);

    string line;
    ps.resize(n);

    for (int i = 0; i < n; i++) {
        double x, y;
        getline(ifs, line);

        istringstream(line) >> x >> y;
        ps[i] = {x, y};
        rs[i][0] = {x, y};
        rs[i][1] = {x, y};
    }

    return true;
}

// 2長方形a,bが交差するか判定 (接するのは許容)
bool is_overlap(const R& a, const R& b) {
    return a[0][0] + E < b[1][0] && b[0][0] + E < a[1][0] &&
           a[0][1] + E < b[1][1] && b[0][1] + E < a[1][1];
}

bool is_overlap(R const& r, P const& p) { return is_overlap(r, {p, p}); }

// 実行可能性のチェック
double area_check(vector<P> const& ps, vector<R> const& rs) {
    // 重なりチェック
    for (int i = 0; i < rs.size() - 1; i++) {
        if (abs(ps[i][0] - rs[i][0][0]) < 10e-6 &&
            abs(ps[i][1] - rs[i][0][1]) < 10e-6) {
            for (int j = i + 1; j < rs.size(); j++) {
                if (is_overlap(rs[i], rs[j])) {
                    // cout << i << " " << j << endl;
                    assert("intersects rs[i] and rs[j]" == "");
                }
            }
        } else {
            assert("anchor is not corR." == "");
        }
    }

    // 面積計算
    double area = 0;
    for (R r : rs) {
        area += (r[1][0] - r[0][0]) * (r[1][1] - r[0][1]);
    }
    return area;
}

array<double, 2> crossover(vector<P> const& ps, vector<R>& rs, mt19937& eng,
                           uniform_real_distribution<double>& dist,
                           vector<V>& gene1, vector<V>& gene2) {
    int n = rs.size();

    vector<vector<V>> children = vector<vector<V>>(CHILDS, vector<V>(n));
    // children[0] = gene1;
    // children[1] = gene2;
    for (int c = 0; c < CHILDS; c++) {
        for (int i = 0; i < n; i++) {
            // 遺伝
            if (dist(eng) < 0.5) {
                children[c][i] = gene1[i];
            } else {
                children[c][i] = gene2[i];
            }

            // 突然変異
            if (dist(eng) < 0.05) {
                children[c][i][0] = dist(eng);
            }

            if (dist(eng) < 0.05) {
                children[c][i][1] = dist(eng);
            }
        }
    }

    // シミュレート
    int max1 = 0, max2 = 0;
    vector<double> result = vector<double>(CHILDS);
    for (int c = 0; c < CHILDS; c++) {
        result[c] = simulate(ps, rs, children[c]);

        if (result[c] > result[max1]) {
            max2 = max1;
            max1 = c;
        } else if (result[c] > result[max2]) {
            max2 = c;
        }
    }

    // 遺伝
    for (int i = 0; i < n; i++) {
        gene1[i] = children[max1][i];
        gene2[i] = children[max2][i];
    }

    return {result[max1], result[max2]};
}

double simulate(vector<P> const& ps, vector<R>& rs, vector<V> const& gene) {
    int n = rs.size();

    // O(n^2)
    // blocks[i][j]: i が j にブロックされる
    // array<bool, 2>:  X方向にブロックされるか, Y方向にブロックされるか
    vector<vector<array<bool, 2>>> blocks = vector<vector<array<bool, 2>>>(
        n, vector<array<bool, 2>>(n, {false, false}));
    vector<vector<double>> times =
        vector<vector<double>>(n, vector<double>(n, -1));
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            P pi = ps[i], pj = ps[j];

            if (pi[0] < pj[0] + E) {
                if (pi[1] < pj[1] + E) {
                    // 第一象限
                    double tx = (pj[0] - pi[0]) / gene[i][0];
                    double ty = (pj[1] - pi[1]) / gene[i][1];
                    if (tx < ty) {
                        // 上がぶつかる
                        blocks[i][j] = {false, true};
                        times[i][j] = ty;
                    } else {
                        // 右がぶつかる
                        blocks[i][j] = {true, false};
                        times[i][j] = tx;
                    }
                    blocks[j][i] = {false, false};
                } else {
                    // 第四象限
                    double ti = (pj[0] - pi[0]) / gene[i][0];
                    double tj = (pi[1] - pj[1]) / gene[j][1];
                    if (ti < tj) {
                        // i優先
                        blocks[i][j] = {false, false};
                        blocks[j][i] = {false, true};
                        times[j][i] = tj;
                    } else {
                        // j優先
                        blocks[i][j] = {true, false};
                        blocks[j][i] = {false, false};
                        times[i][j] = ti;
                    }
                }

            } else {
                if (pi[1] > pj[1] - E) {
                    // 第三象限
                    double tx = (pi[0] - pj[0]) / gene[j][0];
                    double ty = (pi[1] - pj[1]) / gene[j][1];
                    if (tx < ty) {
                        // 上がぶつかる
                        blocks[j][i] = {false, true};
                        times[j][i] = ty;
                    } else {
                        // 右がぶつかる
                        blocks[j][i] = {true, false};
                        times[j][i] = tx;
                    }
                    blocks[i][j] = {false, false};
                } else {
                    // 第二象限
                    double ti = (pj[1] - pi[1]) / gene[i][1];
                    double tj = (pi[0] - pj[0]) / gene[j][0];
                    if (ti < tj) {
                        // i優先
                        blocks[i][j] = {false, false};
                        blocks[j][i] = {true, false};
                        times[j][i] = tj;
                    } else {
                        // j優先
                        blocks[i][j] = {false, true};
                        blocks[j][i] = {false, false};
                        times[i][j] = ti;
                    }
                }
            }
        }
    }

    // 実行結果
    // O(n^2)
    double area = 0;
    bool last_is_right = false;
    for (int i = 0; i < n; i++) {
        double right = 1;
        double top = 1;
        double tr = -1, tt = -1;

        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            if (blocks[i][j][0] && (ps[j][1] - E < top) &&
                (tr < 0 || tr > times[i][j])) {
                tr = times[i][j];
                right = ps[i][0] + gene[i][0] * tr;
                last_is_right = true;
            }

            if (blocks[i][j][1] && (ps[j][0] - E < right) &&
                (tt < 0 || tt > times[i][j])) {
                tt = times[i][j];
                top = ps[i][1] + gene[i][1] * tt;
                last_is_right = false;
            }
        }

        if (last_is_right) {
            // 上に延伸
            double top = 1;
            for (int j = 0; j < n; j++) {
                if (blocks[i][j][1] && (ps[j][0] - E < right) &&
                    (rs[j][1][0] + E > ps[i][0])) {
                    top = min(top, ps[j][1]);
                }
            }
            rs[i][1][1] = top;

        } else {
            // 右に延伸
            double right = 1;
            for (int j = 0; j < n; j++) {
                if (blocks[i][j][0] && (ps[j][1] - E < top) &&
                    (rs[j][1][1] + E > ps[i][1])) {
                    right = min(right, ps[j][0]);
                }
            }
            rs[i][1][0] = right;
        }

        rs[i][1] = {right, top};
        area += rect_area(rs[i]);
    }

    return area;
}