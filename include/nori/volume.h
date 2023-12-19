#if !defined(__NORI_VOLUME_H)
#define __NORI_VOLUME_H
#include <nori/object.h>
#include <nori/warp.h>
#include <nori/medium.h>

NORI_NAMESPACE_BEGIN


class Volume : public NoriObject {
public:
	virtual float lookupFloat(const Point3f& p) const = 0;
	virtual Color3f lookupRGB(const Point3f& p) const = 0;
	EClassType getClassType() const override{ return EVolume; }
	virtual void setMedium(const Medium* medium) { m_medium = medium; }
	virtual bool isPerlin() { return false; }
	virtual void gridGeneration() {}
protected:
	const Medium * m_medium = nullptr;
};

NORI_NAMESPACE_END
#endif /* __NORI_VOLUME_H */
