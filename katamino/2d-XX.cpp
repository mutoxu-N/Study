/*
    カタミノ詰め込みの種類数を計算するプログラム
*/

#include <bits/stdc++.h>

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

using namespace std;

class Coord {
   public:
    int x, y;
};

class Mino {
   public:
    int width, height, type;
    array<Coord, 5> coords;

    bool can_place(vector<vector<char>> const &board, int x, int y) const {
        int board_width = board[0].size(), board_height = board.size();
        if (x + width >= board_width ||
            y + height - coords[0].y >= board_height || y - coords[0].y < 0)
            return false;

        for (int i = 0; i < 5; i++) {
            int nx = x + coords[i].x - coords[0].x;
            int ny = y + coords[i].y - coords[0].y;
            if (board[ny][nx] != -1) return false;
        }
        return true;
    }

    void place(vector<vector<char>> &board, int x, int y) const {
        for (int i = 0; i < 5; i++) {
            int nx = x + coords[i].x - coords[0].x;
            int ny = y + coords[i].y - coords[0].y;
            board[ny][nx] = type;
        }
        return;
    }
};

bool read_minos(vector<vector<Mino>> &minos);
int put(vector<vector<char>> const &board, vector<vector<Mino>> const &remains,
        int start_x, int start_y, int depth, int max_depth);
void print_board(vector<vector<char>> const &board);

vector<vector<vector<char>>> ans;

int main(int argc, char *argv[]) {
    int width = 12, height = 5;
    ans.reserve(10000);

    if (argc > 2) {
        width = atoi(argv[1]);
        height = atoi(argv[2]);
    }

    // if (width < height) {
    //     int t = width;
    //     width = height;
    //     height = t;
    // }

    assert(width * height % 5 == 0);

    vector<vector<Mino>> minos(12);
    assert(read_minos(minos));

    int max_depth = width * height / 5;

    // 時間計測開始
    auto time0 = chrono::system_clock::now();

    // x配置座標
    bool odd_width = (width % 2 == 1);
    bool odd_height = (height % 2 == 1);

    vector<vector<char>> board(height, vector<char>(width, -1));
    int count = put(board, minos, 0, -1, 0, max_depth);

    // 時間計測終了
    auto time1 = chrono::system_clock::now();

    // 整理
    assert(count % 4 == 0);
    count /= 4;

    // 結果出力
    cout << "size: " << width << "x" << height << endl;
    cout << "time: "
         << chrono::duration_cast<chrono::milliseconds>(time1 - time0).count()
         << "ms" << endl;
    cout << "count: " << count << endl;

    return 0;
    // 答え
    int cnt_sub = 0;
    ofstream ofs("ans(2d-XX)_" + to_string(width) + "x" + to_string(height) +
                 ".txt");
    for (int i = 0; i < ans.size(); i++) {
        ofs << i << " / " << ans.size() - 1 << endl;
        vector<vector<char>> const &board = ans[i];
        map<char, string> chars = {
            {F, "F"}, {I, "I"}, {L, "L"},   {N, "N"}, {P, "P"},
            {T, "T"}, {U, "U"}, {V, "V"},   {W, "W"}, {Y, "Y"},
            {Z, "Z"}, {X, "X"}, {EMP, "_"},
        };

        for (vector<char> bs : board) {
            for (char b : bs) {
                ofs << chars[b];
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

void print_board(vector<vector<char>> const &board) {
    map<char, string> chars = {
        {F, "F"}, {I, "I"}, {L, "L"}, {N, "N"}, {P, "P"}, {T, "T"},   {U, "U"},
        {V, "V"}, {W, "W"}, {Y, "Y"}, {Z, "Z"}, {X, "X"}, {EMP, "_"},
    };

    for (vector<char> bs : board) {
        for (char b : bs) {
            cout << chars[b] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int put(vector<vector<char>> const &board, vector<vector<Mino>> const &remains,
        int last_x, int last_y, int depth, int max_depth) {
    if (depth == max_depth) {
        ans.push_back(board);
        return 1;
    }

    int empty_x = -1, empty_y = -1;  // 空きマスの座標

    for (int y = last_y + 1; y < board.size(); y++) {
        if (board[y][last_x] != -1) continue;
        empty_x = last_x;
        empty_y = y;
        break;
    }

    if (empty_x == -1)
        for (int x = last_x + 1; x < board[0].size(); x++) {
            bool f = false;
            for (int y = 0; y < board.size(); y++) {
                if (board[y][x] != -1) continue;
                empty_x = x;
                empty_y = y;
                f = true;
                break;
            }

            if (f) break;
        }
    if (!(empty_x >= 0 && empty_y >= 0)) {
        cout << last_x << " " << last_y << endl;
        cout << empty_x << " " << empty_y << endl;
        cout << depth << " / " << max_depth << endl;
        print_board(board);
    }
    assert(empty_x >= 0 && empty_y >= 0);

    int cnt = 0;
    for (int i = 0; i < remains.size(); i++) {
        vector<Mino> remain = remains[i];
        for (Mino const &mino : remain) {
            if (mino.can_place(board, empty_x, empty_y)) {
                vector<vector<char>> new_board = board;
                for (int j = 0; j < 5; j++) {
                    int nx = empty_x + mino.coords[j].x - mino.coords[0].x;
                    int ny = empty_y + mino.coords[j].y - mino.coords[0].y;
                    new_board[ny][nx] = i;
                }

                vector<vector<Mino>> new_remains = remains;
                new_remains[i].clear();

                cnt += put(new_board, new_remains, empty_x, empty_y, depth + 1,
                           max_depth);
            }
        }
    }
    return cnt;
}
