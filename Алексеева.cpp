#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

using namespace std;

struct Point {
    double x, y;
};

double distance(const Point& a, const Point& b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

double pathLength(const vector<int>& path, const vector<vector<double>>& dist) {
    double length = 0;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        length += dist[path[i]][path[i + 1]];
    }
    length += dist[path.back()][path.front()];
    return length;
}

void two_Opt(vector<int>& path, const vector<vector<double>>& dist) {
    bool improvement = true;
    while (improvement) {
        improvement = false;
        for (size_t i = 1; i < path.size() - 1; ++i) {
            for (size_t j = i + 1; j < path.size(); ++j) {
                if (j == path.size() - 1) continue;
                double delta = -dist[path[i - 1]][path[i]] - dist[path[j]][path[j + 1]]
                             + dist[path[i - 1]][path[j]] + dist[path[i]][path[j + 1]];
                if (delta < 0) {
                    reverse(path.begin() + i, path.begin() + j + 1);
                    improvement = true;
                }
            }
        }
    }
}

void three_Opt(vector<int>& path, const vector<vector<double>>& dist) {
    bool improvement = true;
    while (improvement) {
        improvement = false;
        for (size_t i = 1; i < path.size() - 2; ++i) {
            for (size_t j = i + 1; j < path.size() - 1; ++j) {
                for (size_t k = j + 1; k < path.size(); ++k) {
                    double d1 = dist[path[i - 1]][path[i]] + dist[path[j - 1]][path[j]] + dist[path[k - 1]][path[k]];
                    double d2 = dist[path[i - 1]][path[j]] + dist[path[i]][path[k]] + dist[path[j - 1]][path[k - 1]];
                    if (d2 < d1) {
                        reverse(path.begin() + i, path.begin() + j + 1);
                        reverse(path.begin() + j + 1, path.begin() + k + 1);
                        improvement = true;
                    }
                }
            }
        }
    }
}

int main() {
    ifstream infile("data");
    if (!infile.is_open()) {
        cerr << "Ошибка при открытии файла!" << endl;
        return 1;
    }

    int N;
    infile >> N;
    vector<Point> points(N);
    for (int i = 0; i < N; ++i) {
        infile >> points[i].x >> points[i].y;
    }

    vector<vector<double>> dist(N, vector<double>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            dist[i][j] = distance(points[i], points[j]);
        }
    }

    vector<int> path(N);
    for (int i = 0; i < N; ++i) {
        path[i] = i;
    }

    two_Opt(path, dist);
    three_Opt(path, dist);

    double length = pathLength(path, dist);
    cout << length << " 0" << endl;
    for (int p : path) {
        cout << p << " ";
    }
    cout << endl;

    return 0;
}