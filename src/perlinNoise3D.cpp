#include <iostream>
#include<cstdlib>
#include <cmath>
using namespace std;

class PerlinNoise3D {
public:
    PerlinNoise3D(int row, int col, int height) {
        //generate perlin noise for the grid
        srand((unsigned)time(NULL));
        m_row = row;
        m_col = col;
        m_height = height;
        m_grid = new float[row * col * height];
        for (int k = 0; k < height; k++) {
            for (int j = 0; j < col; j++) {
                for (int i = 0; i < row; i++) {
                    m_grid[i + j * row + k * col * row] = rand();
                }
            }
        }
    }


    /* Given a fixed position, calculate the interpolation value based on its eight neighbour vertex.
       the location is in grid coordinates*/
    virtual float getP(float gridX, float gridY, float gridZ) {
        int beginX = int(floor(gridX));
        int beginY = int(floor(gridY));
        int beginZ = int(floor(gridZ));
        float p1 = interpolate(getValue(beginX, beginY, beginZ), getValue(beginX + 1, beginY, beginZ), lerp(gridX - beginX));
        float p2 = interpolate(getValue(beginX, beginY, beginZ + 1), getValue(beginX + 1, beginY, beginZ + 1), lerp(gridX - beginX));
        float p3 = interpolate(getValue(beginX, beginY + 1, beginZ), getValue(beginX + 1, beginY + 1, beginZ), lerp(gridX - beginX));
        float p4 = interpolate(getValue(beginX, beginY + 1, beginZ + 1), getValue(beginX + 1, beginY + 1, beginZ + 1), lerp(gridX - beginX));
        float p5 = interpolate(p1, p2, lerp(gridZ - beginZ));
        float p6 = interpolate(p3, p4, lerp(gridZ - beginZ));
        return interpolate(p5, p6, lerp(gridY - beginY));
    }

    /* retrieve the grid value*/
    virtual float getValue(int intGridX, int intGridY, int intGridZ) {
        return m_grid[intGridX + intGridY * m_row + intGridZ * m_row * m_col];
    }
    virtual float interpolate(float a, float b, float t) {
        return t * a + (1 - t) * b;
    }

    //perlin interpolation
    virtual float lerp(float num)
    {
        return 6 * pow(num, 5) - 15 * pow(num, 4) + 10 * pow(num, 3);
    }

    virtual float* getgrid() {
        return m_grid;
    }
protected:
    int m_row = 0;
    int m_col = 0;
    int m_height = 0;
    float* m_grid = nullptr;
};