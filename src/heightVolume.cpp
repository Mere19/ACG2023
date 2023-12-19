#include<nori/volume.h>
#include <nori/shape.h>
NORI_NAMESPACE_BEGIN

using namespace std;

/* The Height Volume class procedurally generate the value (primarily density) in the decay of the height*/
class HeightFloatVolume : public Volume {
public:
    HeightFloatVolume(const PropertyList& props) {
        m_value = props.getFloat("value");
    }
    virtual float lookupFloat(const Point3f& p) const override {
        const Shape* mesh = this->m_medium->getShape();
        float zMax = mesh->getZMax();
        float zMin = mesh->getZMin();
        float h = p.z() - zMin;
        float b = 6;
        return m_value * exp(-b * h);
        
    };
    virtual Color3f lookupRGB(const Point3f& p) const override {
        //not used
        return Color3f(0.0f);
    };
    std::string toString() const override {
        return tfm::format("ConstantFloatVolume[]");
    }
protected:
    float m_value;
};

NORI_REGISTER_CLASS(HeightFloatVolume, "heightFVol");
NORI_NAMESPACE_END