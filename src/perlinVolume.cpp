#include<nori/volume.h>
#include <iostream>
#include<cstdlib>
#include <cmath>
#include<nori/shape.h>
NORI_NAMESPACE_BEGIN
using namespace std;

class PerlinNoise3D {
public:
    PerlinNoise3D(int row, int col, int height, float max) {
        srand((unsigned)time(NULL));
        m_row = row;
        m_col = col;
        m_height = height;
        m_grid = new float[row * col * height];
        for (int k = 0; k < height; k++) {
            for (int j = 0; j < col; j++) {
                for (int i = 0; i < row; i++) {
                    m_grid[i + j * row + k * col * row] = static_cast<float>(rand()) / RAND_MAX * max;
                    //cout << m_grid[i + j * row + k * col * row] << endl;
                }
            }
        }
    }

    virtual float getP(float gridX, float gridY, float gridZ) {

        int beginX = int(floor(gridX));
        int beginY = int(floor(gridY));
        int beginZ = int(floor(gridZ));
        //if (beginX < 0 || beginY < 0 || beginZ < 0 || beginX >= m_row - 1 || beginY >= m_col - 1 || beginZ >= m_height - 1) {
        //    std::cerr << "Invalid lookup position or channel" << std::endl;
        //    return 0.0f; // Or any default value as needed
        //}
        float p1 = interpolate(getValue(beginX, beginY, beginZ), getValue(beginX + 1, beginY, beginZ), lerp(gridX - beginX));
        float p2 = interpolate(getValue(beginX, beginY, beginZ + 1), getValue(beginX + 1, beginY, beginZ + 1), lerp(gridX - beginX));
        float p3 = interpolate(getValue(beginX, beginY + 1, beginZ), getValue(beginX + 1, beginY + 1, beginZ), lerp(gridX - beginX));
        float p4 = interpolate(getValue(beginX, beginY + 1, beginZ + 1), getValue(beginX + 1, beginY + 1, beginZ + 1), lerp(gridX - beginX));
        float p5 = interpolate(p1, p2, lerp(gridZ - beginZ));
        float p6 = interpolate(p3, p4, lerp(gridZ - beginZ));
        return interpolate(p5, p6, lerp(gridY - beginY));
    }
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

class PerlinVolume : public Volume {
public:
    PerlinVolume(const PropertyList& props) {
        m_value = props.getFloat("value");
    }
    virtual void gridGeneration() override {
        
        BoundingBox3f bbox = this->m_medium->getShape()->getBoundingBox();
        this->m_medium->getShape();

        m_min = bbox.min;
        m_max = bbox.max;
        int row = (int)ceil((m_max.x() - m_min.x()) / m_scale_x);
        float col = (int)ceil((m_max.y() - m_min.y()) / m_scale_y);
        float height = (int)ceil((m_max.z() - m_min.z()) / m_scale_z);
        m_pn = new PerlinNoise3D(row, col, height, m_value);
    }
    virtual float lookupFloat(const Point3f& p) const override {
        float result = m_pn->getP((p.x()-m_min.x())/m_scale_x, (p.y() - m_min.y()) / m_scale_y, (p.z() - m_min.z()) / m_scale_z);
        return result;
    };
    virtual Color3f lookupRGB(const Point3f& p) const override {
        //not used
        return Color3f(0.0f);
    };
    virtual bool isPerlin() override{
        return true;
    }
    std::string toString() const override {
        return tfm::format("PerlinFloatVolume[]");
    }
protected:
    float m_value;
    float m_scale_x = 0.1;
    float m_scale_y = 0.1;
    float m_scale_z = 0.1;
    Point3f m_min = Point3f(0);
    Point3f m_max = Point3f(0);
    PerlinNoise3D *m_pn = nullptr;
};

NORI_REGISTER_CLASS(PerlinVolume, "plVol");
NORI_NAMESPACE_END