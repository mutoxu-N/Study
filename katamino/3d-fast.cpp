/*
    カタミノ3D詰め込みの種類数を計算するプログラム
    - Iミノの位置を考慮することで効率化
    - OpenMPを使った並列化
*/

#include <bits/stdc++.h>
#include <omp.h>

#define EMP -1
#define F 0
#define I 1
#define L 2
#define N 3
#define P 4
#define T 5
#define U 6
#define V 7
#define W 8
#define Y 9
#define Z 10
#define X 11
#define D1 5
#define D2 4
#define D3 3

using namespace std;

class Coord {
   public:
    int x, y;
};

class Mino {
   public:
    int width, height, type;
    array<Coord, 5> coords;

    bool can_place(vector<vector<vector<char>>> const &board, int x, int y,
                   int z, int d) const {
        switch (d) {
            case 0:
                // xy平面
                if (x + width >= D1 || y + height - coords[0].y >= D2 ||
                    y - coords[0].y < 0)
                    return false;

                for (int i = 0; i < 5; i++) {
                    int nx = x + coords[i].x - coords[0].x;
                    int ny = y + coords[i].y - coords[0].y;
                    if (board[nx][ny][z] != -1) return false;
                }
                break;
            case 1:
                // yz平面
                if (y + width >= D2 || z + height - coords[0].y >= D3 ||
                    z - coords[0].y < 0)
                    return false;

                for (int i = 0; i < 5; i++) {
                    int ny = y + coords[i].x - coords[0].x;
                    int nz = z + coords[i].y - coords[0].y;
                    if (board[x][ny][nz] != -1) return false;
                }
                break;
            case 2:
                // xz平面
                if (x + width >= D1 ||
                    z + height - coords[0].y >= z - coords[0].y < 0)
                    return false;

                for (int i = 0; i < 5; i++) {
                    int nx = x + coords[i].x - coords[0].x;
                    int nz = z + coords[i].y - coords[0].y;
                    if (board[nx][y][nz] != -1) return false;
                }
                break;
        }

        return true;
    }

    void place(vector<vector<vector<char>>> &board, int x, int y, int z,
               int d) const {
        switch (d) {
            case 0:
                // xy平面
                for (int i = 0; i < 5; i++) {
                    int nx = x + coords[i].x - coords[0].x;
                    int ny = y + coords[i].y - coords[0].y;
                    board[nx][ny][z] = type;
                }
                break;
            case 1:
                // yz平面
                for (int i = 0; i < 5; i++) {
                    int ny = y + coords[i].x - coords[0].x;
                    int nz = z + coords[i].y - coords[0].y;
                    board[x][ny][nz] = type;
                }
                break;
            case 2:
                // xz平面
                for (int i = 0; i < 5; i++) {
                    int nx = x + coords[i].x - coords[0].x;
                    int nz = z + coords[i].y - coords[0].y;
                    board[nx][y][nz] = type;
                }
                break;
        }
    }
};

bool read_minos(vector<vector<Mino>> &minos);
int put(vector<vector<vector<char>>> const &board,
        vector<vector<Mino>> const &remains, int last, int depth,
        int max_depth);
void print_board(vector<vector<vector<char>>> const &board);

vector<vector<vector<vector<char>>>> ans;

int main(int argc, char *argv[]) {
    ans.reserve(10000);

    vector<vector<Mino>> minos(12);
    assert(read_minos(minos));
    minos[I].clear();

    int max_depth = 12;

    // 時間計測開始
    auto time0 = chrono::system_clock::now();

    // B: x2
    array<int, 2> counts = {0, 0};
    for (int y = 0; y < 2; y++) {
#pragma omp parallel for
        for (int x = 0; x < 2; x++) {
            vector<vector<vector<char>>> board(
                D1, vector<vector<char>>(D2, vector<char>(D3, -1)));

            for (int z = 0; z < 5; z++) {
                board[z][x][y] = I;
            }
#pragma omp critical
            counts[y] += put(board, minos, 0, 1, max_depth);
        }
    }

    // 時間計測終了
    auto time1 = chrono::system_clock::now();

    // 結果整理
    int count_a = counts[0], count_b = counts[1];
    int count = count_a / 2 + count_b / 4;
    // assert(count_a % 2 == 0);
    // assert(count_b % 4 == 0);

    // 結果出力
    cout << "size: " << D1 << "x" << D2 << "x" << D3 << endl;
    cout << "time: "
         << chrono::duration_cast<chrono::milliseconds>(time1 - time0).count()
         << "ms" << endl;
    cout << "count: " << count << endl;
    cout << "\ta: " << count_a << endl;
    cout << "\tb: " << count_b << endl;

    // 答え
    int cnt_sub = 0;
    ofstream ofs("ans(3d-fast).txt");
    for (int i = 0; i < ans.size(); i++) {
        ofs << i << " / " << ans.size() - 1 << endl;
        vector<vector<vector<char>>> const &board = ans[i];
        map<char, string> chars = {
            {F, "F"}, {I, "I"}, {L, "L"},   {N, "N"}, {P, "P"},
            {T, "T"}, {U, "U"}, {V, "V"},   {W, "W"}, {Y, "Y"},
            {Z, "Z"}, {X, "X"}, {EMP, "_"},
        };

        for (int z = 0; z < D3; z++) {
            for (int y = 0; y < D2; y++) {
                for (int x = 0; x < D1; x++) {
                    ofs << chars[board[x][y][z]];
                }
                ofs << endl;
            }
            ofs << endl;
        }
        ofs << endl;
    }

    return 0;
}

bool read_minos(vector<vector<Mino>> &minos) {
    ifstream ifs("minos.txt");

    if (!ifs) {
        cerr << "file \"minos.txt\" cannot be found!" << endl;
        return false;
    }

    int w, h, x, y;
    string line, t;

    while (getline(ifs, line)) {
        if (line.empty()) continue;

        istringstream iss(line);
        iss >> t >> w >> h;

        array<Coord, 5> coords;

        for (int i = 0; i < 5; i++) {
            getline(ifs, line);
            istringstream iss(line);

            iss >> x >> y;
            coords[i] = Coord{x, y};
        }

        char idx = -1;
        switch (t[0]) {
            case 'x':
                idx += 1;
            case 'z':
                idx += 1;
            case 'y':
                idx += 1;
            case 'w':
                idx += 1;
            case 'v':
                idx += 1;
            case 'u':
                idx += 1;
            case 't':
                idx += 1;
            case 'p':
                idx += 1;
            case 'n':
                idx += 1;
            case 'l':
                idx += 1;
            case 'i':
                idx += 1;
            case 'f':
                idx += 1;
        }

        assert(idx != -1);
        minos[idx].push_back(Mino{w, h, idx, coords});
    }

    return true;
}

void print_board(vector<vector<vector<char>>> const &board) {
    map<char, string> chars = {
        {F, "F"}, {I, "I"}, {L, "L"}, {N, "N"}, {P, "P"}, {T, "T"},   {U, "U"},
        {V, "V"}, {W, "W"}, {Y, "Y"}, {Z, "Z"}, {X, "X"}, {EMP, "_"},
    };

    for (int z = 0; z < board[0][0].size(); z++) {
        for (int y = 0; y < board[0].size(); y++) {
            for (int x = 0; x < board.size(); x++) {
                cout << chars[board[x][y][z]] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}

int put(vector<vector<vector<char>>> const &board,
        vector<vector<Mino>> const &remains, int last, int depth,
        int max_depth) {
    if (depth == max_depth) {
        ans.push_back(board);
        // cout << "found!: " << ans.size() << endl;
        return 1;
    }

    int empty_x = -1, empty_y = -1, empty_z = -1;
    for (int i = last; i < 60; i++) {
        int x = i / (D2 * D3);
        int y = i % (D2 * D3) / D3;
        int z = i % D3;

        if (board[x][y][z] != EMP) continue;

        empty_x = x;
        empty_y = y;
        empty_z = z;
        break;
    }

    if (!(empty_x != -1 && empty_y != -1 && empty_z != -1)) {
        cout << "last: " << last << endl;
        cout << "depth: " << depth << endl;
        print_board(board);
    }
    assert(empty_x != -1 && empty_y != -1 && empty_z != -1);

    int cnt = 0;
    for (int i = 0; i < remains.size(); i++) {
        vector<Mino> remain = remains[i];
        for (Mino const &mino : remain) {
            for (int d = 0; d < 3; d++) {
                if (mino.can_place(board, empty_x, empty_y, empty_z, d)) {
                    vector<vector<vector<char>>> new_board = board;
                    mino.place(new_board, empty_x, empty_y, empty_z, d);

                    vector<vector<Mino>> new_remains = remains;
                    new_remains[i].clear();

                    cnt += put(new_board, new_remains,
                               (empty_x * D3 * D2 + empty_y * D3 + empty_z),
                               depth + 1, max_depth);
                }
            }
        }
    }
    return cnt;
}
