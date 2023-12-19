#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <nori/volume.h>
#include <nori/shape.h>

NORI_NAMESPACE_BEGIN
using namespace std;
/* The Grid Volume class parse a Mitsuba-standard .vol file, create a corresponding axis-aligned boundingbox.
and is responsible for retrieving the value at a given point based on its position inside the bounding box*/

class GridVolume : public Volume {

public:
    GridVolume(const PropertyList& props) {
        m_filename =props.getString("filename");
        std::ifstream file(m_filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << m_filename << std::endl;
            return;
        }

        //byte 1-3
        char header[3];
        file.read(header, 3);
        if (std::string(header, 3) != "VOL") {
            std::cerr << "Invalid volume data file!" << std::endl;
            return;
        }

        //byte 4
        int version;
        file.read(reinterpret_cast<char*>(&version), sizeof(char));
        //if (version != 3) {
        //    std::cerr << "Unsupported version!" << std::endl;
        //    return;
        //}
        //byte 5-8
        int encoding;
        file.read(reinterpret_cast<char*>(&encoding), sizeof(int));
        //byte 9-12
        file.read(reinterpret_cast<char*>(&m_xres), sizeof(int));
        //byte 13-16
        file.read(reinterpret_cast<char*>(&m_yres), sizeof(int));
        //byte 17-20
        file.read(reinterpret_cast<char*>(&m_zres), sizeof(int));
        //byte 21-24
        file.read(reinterpret_cast<char*>(&m_channels), sizeof(int));
        cout << "channel number: " << m_channels << endl;
        cout << "xres: " << m_xres << endl;
        cout << "yres: " << m_yres << endl;
        cout << "zres: " << m_zres << endl;
        //byte 25-48
        float boundingBox[6];
        for (int i = 0; i < 6; ++i) {
            file.read(reinterpret_cast<char*>(&boundingBox[i]), sizeof(float));
            cout << "bb" << i << ": " << boundingBox[i] << endl;
        }
        //Read binary data
        int dataSize = m_xres * m_yres * m_zres * m_channels;
        data.resize(dataSize);
        file.read(reinterpret_cast<char*>(data.data()), dataSize * sizeof(float));
        file.close();
    }

    float lookupFloat(const Point3f& p) const override{
        BoundingBox3f bbox = this->m_medium->getShape()->getBoundingBox();
        Vector3f diag = bbox.getExtents();
        Point3f min = bbox.min;
        Point3f bb_p = p - min;
        bb_p = Point3f(bb_p.x() / diag.x() * (m_xres - 1), bb_p.y() / diag.y() * (m_yres - 1), bb_p.z() / diag.z() * (m_zres - 1));
        int x1 = static_cast<int>(std::floor(bb_p.x()));
        int y1 = static_cast<int>(std::floor(bb_p.y()));
        int z1 = static_cast<int>(std::floor(bb_p.z()));
        int chan = 0;
        if (x1 < 0 || y1 < 0 || z1 < 0 || x1 >= m_xres - 1 || y1 >= m_yres - 1 || z1 >= m_zres - 1 || chan < 0 || chan >= m_channels) {
            std::cerr << "Invalid lookup position or channel" << std::endl;
            return 0.0f; // Or any default value as needed
        }

        float xd = bb_p.x() - x1;
        float yd = bb_p.y() - y1;
        float zd = bb_p.z() - z1;


        //trilinear interpolation to get the value
        float c000 = data[((z1 * m_yres + y1) * m_xres + x1) * m_channels + chan];
        float c001 = data[((z1 * m_yres + y1) * m_xres + x1 + 1) * m_channels + chan];
        float c010 = data[((z1 * m_yres + y1 + 1) * m_xres + x1) * m_channels + chan];
        float c011 = data[((z1 * m_yres + y1 + 1) * m_xres + x1 + 1) * m_channels + chan];
        float c100 = data[((z1 + 1) * m_yres + y1) * m_xres + x1 * m_channels + chan];
        float c101 = data[((z1 + 1) * m_yres + y1) * m_xres + (x1 + 1) * m_channels + chan];
        float c110 = data[((z1 + 1) * m_yres + y1 + 1) * m_xres + x1 * m_channels + chan];
        float c111 = data[((z1 + 1) * m_yres + y1 + 1) * m_xres + (x1 + 1) * m_channels + chan];

        float c00 = c000 * (1 - xd) + c001 * xd;
        float c01 = c010 * (1 - xd) + c011 * xd;
        float c10 = c100 * (1 - xd) + c101 * xd;
        float c11 = c110 * (1 - xd) + c111 * xd;

        float c0 = c00 * (1 - yd) + c01 * yd;
        float c1 = c10 * (1 - yd) + c11 * yd;

        return c0 * (1 - zd) + c1 * zd;
    }
    virtual Color3f lookupRGB(const Point3f& p) const override {
        Color3f result(0);
        BoundingBox3f bbox = this->m_medium->getShape()->getBoundingBox();
        Vector3f diag = bbox.getExtents();
        Point3f min = bbox.min;
        Point3f bb_p = p - min;
        bb_p = Point3f(bb_p.x() / diag.x() * (m_xres - 1), bb_p.y() / diag.y() * (m_yres - 1), bb_p.z() / diag.z() * (m_zres - 1));
        int x1 = static_cast<int>(std::floor(bb_p.x()));
        int y1 = static_cast<int>(std::floor(bb_p.y()));
        int z1 = static_cast<int>(std::floor(bb_p.z()));
        for (int i = 0; i < m_channels; i++) {
            if (x1 < 0 || y1 < 0 || z1 < 0 || x1 >= m_xres - 1 || y1 >= m_yres - 1 || z1 >= m_zres - 1 || i < 0 || i >= m_channels) {
                std::cerr << "Invalid lookup position or channel" << std::endl;
                return 0.0f; // Or any default value as needed
            }

            float xd = bb_p.x() - x1;
            float yd = bb_p.y() - y1;
            float zd = bb_p.z() - z1;

            float c000 = data[((z1 * m_yres + y1) * m_xres + x1) * m_channels + i];
            float c001 = data[((z1 * m_yres + y1) * m_xres + x1 + 1) * m_channels + i];
            float c010 = data[((z1 * m_yres + y1 + 1) * m_xres + x1) * m_channels + i];
            float c011 = data[((z1 * m_yres + y1 + 1) * m_xres + x1 + 1) * m_channels + i];
            float c100 = data[((z1 + 1) * m_yres + y1) * m_xres + x1 * m_channels + i];
            float c101 = data[((z1 + 1) * m_yres + y1) * m_xres + (x1 + 1) * m_channels + i];
            float c110 = data[((z1 + 1) * m_yres + y1 + 1) * m_xres + x1 * m_channels + i];
            float c111 = data[((z1 + 1) * m_yres + y1 + 1) * m_xres + (x1 + 1) * m_channels + i];

            float c00 = c000 * (1 - xd) + c001 * xd;
            float c01 = c010 * (1 - xd) + c011 * xd;
            float c10 = c100 * (1 - xd) + c101 * xd;
            float c11 = c110 * (1 - xd) + c111 * xd;

            float c0 = c00 * (1 - yd) + c01 * yd;
            float c1 = c10 * (1 - yd) + c11 * yd;

            result[i] = c0 * (1 - zd) + c1 * zd;
        }
        return result;
    }
    std::string toString() const override {
        return tfm::format("GridVolume[]");
    }
private:
    std::vector<float> data;
    int m_xres, m_yres, m_zres, m_channels;
    std::string m_filename;
};
NORI_REGISTER_CLASS(GridVolume, "grid");
NORI_NAMESPACE_END