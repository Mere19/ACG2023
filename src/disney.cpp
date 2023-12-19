#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>
#include <nori/common.h>

#include <Eigen/Geometry>
#include <Eigen/LU>

NORI_NAMESPACE_BEGIN

Vector3f cross(const Vector3f &a, const Vector3f &b)
{
    Vector3f r(0.f);
    r(0) = a(1)*b(2)-a(2)*b(1);
    r(1) = a(2)*b(0)-a(0)*b(2);
    r(2) = a(0)*b(1)-a(1)*b(0);
    return r;
}

class Disney : public BSDF {
public:
    Disney(const PropertyList &propList) : m_albedo(nullptr) {
        if(propList.has("albedo")) {
            PropertyList l;
            l.setColor("value", propList.getColor("albedo"));
            m_albedo = static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l));
        }
        cout << m_albedo << endl;

        m_subsurface = propList.getFloat("subsurface", 0.f);
        m_sheen = propList.getFloat("sheen", 0.f);
        m_sheenTint = propList.getFloat("sheenTint", 0.f);
        m_clearcoat = propList.getFloat("clearcoat", 0.f);
        m_clearcoatGloss = propList.getFloat("clearcoatGloss", 0.f);
        m_specular = propList.getFloat("specular", 0.f);
        m_specTint = propList.getFloat("specTint", 0.f);
        m_specTrans = propList.getFloat("specTrans", 0.f);
        m_roughness = propList.getFloat("roughness", 0.f);
        m_anisotropic = propList.getFloat("anistropic", 0.f);
        m_metallic = propList.getFloat("metallic", 0.f);
        m_eta = 2.f / (1.f - sqrt(0.08f * m_specular)) - 1.f;

        // Secondary parameters
        m_alpha_min = 1e-4;
        m_aspect = sqrt(1.f - 0.9 * m_anisotropic);
        m_alphax = fmax(0.0001, m_roughness * m_roughness / m_aspect);
        m_alphay = fmax(0.0001, m_roughness * m_roughness * m_aspect);
        m_alpha_g = (1 - m_clearcoatGloss) * 0.1f + m_clearcoatGloss * 0.001;
    }
    virtual ~Disney() {
        delete m_albedo;
    }

    /// Add texture for the albedo
    virtual void addChild(NoriObject *obj) override {
        switch (obj->getClassType()) {
            case ETexture:
                if(obj->getIdName() == "albedo") {
                    if (m_albedo)
                        throw NoriException("There is already an albedo defined!");
                    m_albedo = static_cast<Texture<Color3f> *>(obj);
                }
                else {
                    throw NoriException("The name of this texture does not match any field!");
                }
                break;

            default:
                throw NoriException("Diffuse::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }

    virtual void activate() override {
        if(!m_albedo) {
            PropertyList l;
            l.setColor("value", Color3f(0.5f));
            m_albedo = static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l));
            m_albedo->activate();
        }
    }

    /* Diffuse 
    *  Parameters:
    *   baseColor
    *   roughness
    *   subsurface
    */
    Color3f evalBase(const BSDFQueryRecord &bRec) const {
        /* compute FD90 */
        float FD90 = 0.5f + 2 * m_roughness * abs(bRec.wo.dot((bRec.wi + bRec.wo).normalized()));

        /* compute FD */
        float FD_wi = 1 + (FD90 - 1) * (1 - pow(abs(bRec.wi.z()), 5));

        return m_albedo->eval(bRec.uv) * INV_PI * FD_wi * FD90 * abs(bRec.wo.z());
    }

    Color3f evalSubsurface(const BSDFQueryRecord &bRec) const {
        /* compute FSS90 */
        float FSS90 = m_roughness * pow(abs(bRec.wo.dot((bRec.wi + bRec.wo).normalized())), 2);

        /* compute FSS_wi */
        float FSS_wi = (1 + (FSS90 - 1) * pow(1 - abs(bRec.wi.z()), 5));

        /* compute FSS_wo */
        float FSS_wo = (1 + (FSS90 - 1) * pow(1 - abs(bRec.wo.z()), 5));

        Color3f x = 1.25f * m_albedo->eval(bRec.uv) * INV_PI;
        Color3f y = FSS_wi * FSS_wo * (1.f / (abs(bRec.wi.z()) + abs(bRec.wo.z())) - 0.5f) + 0.5f;
        Color3f z = abs(bRec.wo.z());

        return x * y * z;
    }

    Color3f evalDiffuse(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        Color3f base = evalBase(bRec);
        Color3f subsurface = evalSubsurface(bRec);

        return (1.f - m_subsurface) * base + m_subsurface * subsurface;
    }

    float pdfDiffuse(const BSDFQueryRecord &bRec) const
    {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;
        
        return INV_PI * Frame::cosTheta(bRec.wo);
    }

    Color3f sampleDiffuse(BSDFQueryRecord &bRec, const Point2f &sample) const {
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;

        /* Warp a uniformly distributed sample on [0,1]^2
           to a direction on a cosine-weighted hemisphere */
        bRec.wo = Warp::squareToCosineHemisphere(sample);

        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;

        /* eval() / pdf() * cos(theta) = albedo. There
           is no need to call these functions. */
        return eval(bRec) / pdf(bRec) * Frame::cosTheta(bRec.wo);
    }


    /* Metallic
    *  Parameters:
    *   baseColor
    *   roughness
    *   anisotropic
    */
    Color3f Fresnel(const BSDFQueryRecord& bRec, const Vector3f& wh) const {
        Color3f basecolor = m_albedo->eval(bRec.uv);
        float luminance = basecolor.getLuminance();
        Color3f Ctint = luminance > 0 ? basecolor / luminance : Color3f(1.f);
        float R0_eta = (m_eta - 1) * (m_eta - 1) / (m_eta + 1) / (m_eta + 1);
        Color3f Ks = (1 - m_specTint) + m_specTint * Ctint;
        Color3f C0 = m_specular * R0_eta * (1 - m_metallic) * Ks + m_metallic * basecolor;

        return C0 + (1.f - C0) * float(pow(1 - wh.dot(bRec.wo), 5));
    }

    /// Evaluate the Trowbridge-Reitz distribution Dm
    /// ref: pbrb - microfacet
    float TrowbridgeReitz(const Vector3f& wh) const {
        // parallel, no sampling
        if (Frame::cosTheta(wh) == 0.f) {
            return 0.f;
        }

        float tan2Theta = Frame::tanTheta(wh) * Frame::tanTheta(wh);
        float cos2Theta = Frame::cosTheta(wh) * Frame::cosTheta(wh);
        float cos4Theta = cos2Theta * cos2Theta;

        float e = (Frame::cosPhi2(wh) / (m_alphax * m_alphax) +
                    Frame::sinPhi2(wh) / (m_alphay * m_alphay)) * tan2Theta;
        
        return INV_PI / (m_alphax * m_alphay * cos4Theta * (1 + e) * (1 + e));
    }

    /// Evaluate the Lambda function
    /// ref: pbrb - microfacet
    float Lambda(const Vector3f& w) const {
        // parallel, no microfacet seen
        if (Frame::cosTheta(w) == 0.f) {
            return 0.f;
        }

        float absTanTheta = abs(Frame::tanTheta(w));
        float cosPhi = Frame::cosPhi(w);
        float sinPhi = Frame::sinPhi(w);

        float alpha2 = cosPhi * cosPhi * m_alphax * m_alphax + sinPhi * sinPhi * m_alphay * m_alphay;
        float alpha2Tan2Theta = alpha2 * absTanTheta * absTanTheta;

        return (-1.f + sqrt(1.f + alpha2Tan2Theta)) / 2.f;
    }

    /// Evaluate Smith's shadowing-masking function G
    float SmithG(const Vector3f &w) const {
        return 1.f / (1.f + Lambda(w));
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f evalMetallic(const BSDFQueryRecord &bRec) const {
        /* no reflection from the back */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        Color3f Fm = Fresnel(bRec, wh);
        float Dm = TrowbridgeReitz(wh);
        float Gm = SmithG(bRec.wi) * SmithG(bRec.wo);
        float cosTheta = Frame::cosTheta(bRec.wi);

        return Fm * Dm * Gm / (4 * cosTheta);
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    /// ref: https://cseweb.ucsd.edu/~tzli/cse272/wi2023/homework1.pdf
    float pdfMetallic(const BSDFQueryRecord &bRec) const {
        if (Frame::cosTheta(bRec.wo) <= 0.f || bRec.measure != ESolidAngle)
            return 0;

        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        float D = TrowbridgeReitz(wh);
        float Gi = SmithG(bRec.wi);
        float Jh = 1.f / (4.f * abs(Frame::cosTheta(bRec.wi)));

        return D * Gi * Jh;
    }

    /// VNDF sample the BRDF
    /// ref: https://jcgt.org/published/0007/04/01/
    Color3f sampleMetallic(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        bRec.eta = 1.f;
        bRec.measure = ESolidAngle;
        
        if (Frame::cosTheta(bRec.wi) <= 0.f)
            return 0;

        // transform the view direction to the hemisphere configuration
        Vector3f Vh = Vector3f(m_alphax * bRec.wi.x(), m_alphay * bRec.wi.y(), bRec.wi.z()).normalized();

        // orthonormal basis
        float lensq = Vh.x() * Vh.x() + Vh.y() * Vh.y();
        Vector3f T1 = lensq > 0 ? Vector3f(-Vh.y(), Vh.x(), 0) / sqrt(lensq) : Vector3f(1, 0, 0);
        Vector3f T2 = Vh.cross(T1);

        // parameterization
        float r = sqrt(_sample.x());
        float phi = 2.f * M_PI * _sample.y();
        float t1 = r * cos(phi);
        float t2 = r * sin(phi);
        float s = 0.5 * (1.0 + Vh.z());
        t2 = (1.0 - s) * sqrt(1.0 - t1 * t1) + s * t2;

        // reprojection onto hemisphere
        Vector3f Nh = t1 * T1 + t2 * T2 + sqrt(fmax(0.0, 1.0 - t1*t1 - t2*t2)) * Vh;

        // transforming normal back to ellipsoid configuration
        Vector3f Ne = Vector3f(m_alphax * Nh.x(), m_alphay * Nh.y(), fmax(0.0, Nh.z())).normalized();
        bRec.wo = ((2.f * bRec.wi.dot(Ne) * Ne) - bRec.wi).normalized();

        /* discard invalid samples */
        if (pdf(bRec) > 0) {
            return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);   // cos forshortening term
        } else {
            return 0;
        }
    }

    /* Clearcoat 
    *  Parameters:
    *    clearcoatGloss
    */
   float R_0(float eta) const
    {
        return (eta - 1.f) * (eta - 1.f) / (eta + 1.f) / (eta + 1.f);
    }

    float Lambda_c_w(Vector3f w) const
    {
        return ( sqrt( 1.f + (pow(w.x() * 0.25f, 2) + pow(w.y() * 0.25f, 2) ) / pow(w.z(), 2) ) - 1.f ) / 2.f;
    }

    float G_c_w(Vector3f w) const
    {
        return 1.f / (1 + Lambda_c_w(w));
    }

    Color3f evalClearcoat(const BSDFQueryRecord &bRec) const
    {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);
            
        Vector3f wh = (bRec.wi + bRec.wo).normalized();

        float x = clamp(1 - abs(wh.dot(bRec.wo)), 0.f, 1.f);
        float Fc = R_0(1.5) + (1 - R_0(1.5)) * pow(x, 5);
        
        float tmp = M_PI * log(m_alpha_g * m_alpha_g) * (1 + (m_alpha_g * m_alpha_g -1) * wh.z() * wh.z());
        float Dc = (m_alpha_g * m_alpha_g - 1) / tmp;

        float Gc = G_c_w(bRec.wi) * G_c_w(bRec.wo);

        return Fc * Dc * Gc / (4 * abs(Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo)));
    }

    float pdfClearcoat(const BSDFQueryRecord &bRec) const
    {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        float tmp = M_PI * log(m_alpha_g * m_alpha_g) * (1 + (m_alpha_g * m_alpha_g -1) * wh.z() * wh.z());
        float Dc = (m_alpha_g * m_alpha_g - 1) / tmp;

        return Dc * abs(Frame::cosTheta(wh));
    }

    void sampleClearcoat(BSDFQueryRecord &bRec, const Point2f &sample) const
    {
        float h_elevation = sqrt( (1 - pow(m_alpha_g * m_alpha_g, 2), 1 - sample.x()) / (1 - m_alpha_g * m_alpha_g));
        h_elevation = acos(clamp(h_elevation, -1.f, 1.f));
        float h_azimuth = 2 * M_PI * sample.y();
        Vector3f wh;
        wh.x() = sin(h_elevation) * cos(h_azimuth);
        wh.y() = sin(h_elevation) * sin(h_azimuth);
        wh.z() = cos(h_elevation);
        wh.normalize();
        bRec.wo = 2.f*(bRec.wi.dot(wh)) * wh - bRec.wi;
    }

    /* Sheen */
    /// Evaluate the Fresnel term Fm
    Color3f evalSheen(const BSDFQueryRecord& bRec) const {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);
            
        Color3f basecolor = m_albedo->eval(bRec.uv);
        float luminance = basecolor.getLuminance();

        Color3f Ctint = luminance > 0 ? basecolor / luminance : Color3f(1.f);
        Color3f Csheen = (1 - m_sheenTint) + m_sheenTint * Ctint;
        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        
        return Csheen * pow((1.f - abs(wh.dot(bRec.wo))), 5) * Frame::cosTheta(bRec.wo);
    }

    float pdfSheen(const BSDFQueryRecord &bRec) const
    {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;
        
        return INV_PI * Frame::cosTheta(bRec.wo);
    }

    /// Sample the BRDF
    virtual Color3f sampleSheen(BSDFQueryRecord &bRec, const Point2f &sample) const {
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;

        /* Warp a uniformly distributed sample on [0,1]^2
           to a direction on a cosine-weighted hemisphere */
        bRec.wo = Warp::squareToCosineHemisphere(sample);

        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;

        /* eval() / pdf() * cos(theta) = albedo. There
           is no need to call these functions. */
        return eval(bRec) / pdf(bRec) * Frame::cosTheta(bRec.wo);
    }

    /// Evaluate the BRDF model
    virtual Color3f eval(const BSDFQueryRecord &bRec) const override {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        Color3f principledBSDF(0.f);

        
        principledBSDF =  (1.f - m_metallic) * evalDiffuse(bRec) 
                        + (1.f - m_metallic) * m_sheen * evalSheen(bRec) 
                        + evalMetallic(bRec)
                        + 0.25 * m_clearcoat * evalClearcoat(bRec);

        return principledBSDF;
    }

    /// Compute the density of \ref sample() wrt. solid angles
    virtual float pdf(const BSDFQueryRecord &bRec) const override {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        float principledPDF;

        float diffuseWeight = 1.f - m_metallic;
        float metalWeight = 1.f;
        float clearcoatWeight = 0.25 * m_clearcoat;

        Vector3f pdf(diffuseWeight, metalWeight, clearcoatWeight);
        pdf.normalize();
        
        Vector3f wh = (bRec.wi, bRec.wo).normalized();
        principledPDF =   pdf(0) * pdfDiffuse(bRec) 
                        + pdf(1) * pdfMetallic(bRec) / (4.f * bRec.wo.dot(wh));
                        + pdf(2) * pdfClearcoat(bRec) / (4.f * bRec.wo.dot(wh));

        return principledPDF;
    }

    /// Draw a sample from the BRDF model
    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const override {
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;

        float diffuseWeight = (1.f - m_metallic * (1 - m_specTrans));
        float metalWeight =  (1 - m_specTrans * (1 - m_metallic));
        float clearcoatWeight = 0.25 * m_clearcoat;

        // Construct cdf
        Vector3f cdf(diffuseWeight, metalWeight, clearcoatWeight);
        cdf /= cdf.sum();
        cdf.y() = cdf.x() + cdf.y();
        cdf.z() = cdf.y() + cdf.z();
        // cout << cdf << endl;

        Point2f sample = _sample;

        if (sample.x() < cdf.x()) // Diffuse
        {
            sample.x() /=  cdf.x();
            sampleDiffuse(bRec, sample);
        }
        else if(sample.x() < cdf.y()) // Metallic
        {
            sample.x() = (sample.x() - cdf.x()) / (cdf.y() - cdf.x());
            sampleMetallic(bRec, sample);
        }
        else // Clearcoat
        {
            sample.x() = (sample.x() - cdf.y()) / (1 - cdf.y());
            sampleClearcoat(bRec, sample);
        }

        float pdf_disney = pdf(bRec);

        if (pdf_disney > 0)
            return eval(bRec) / pdf_disney * std::max(Frame::cosTheta(bRec.wo), 0.f);
        else
            return Color3f(0.f);
    }

    virtual bool isDiffuse() const override {return true;}

    /// Return a human-readable summary
    virtual std::string toString() const override {
        return tfm::format(
            "Disney[\n"
            "  baseColor = %s,\n"
            "  subsurface = %f,\n"
            "  sheen = %f,\n"
            "  sheenTint = %f,\n"
            "  clearcoat = %f,\n"
            "  clearcoatGloss = %f,\n"
            "  specular = %f,\n"
            "  specTint= %f,\n"
            "  specTrans= %f,\n"
            "  roughness = %f,\n"
            "  anistropic = %f,\n"
            "  metallic = %f\n"
            "]",
            m_albedo ? indent(m_albedo->toString()) : std::string("null"),
            m_subsurface,
            m_sheen,
            m_sheenTint,
            m_clearcoat,
            m_clearcoatGloss,
            m_specular,
            m_specTint,
            m_specTrans,
            m_roughness,
            m_anisotropic,
            m_metallic
        );
    }

    virtual EClassType getClassType() const override { return EBSDF; }

private:
    Texture<Color3f> * m_albedo;
    float m_subsurface;
    float m_sheen;
    float m_sheenTint;
    float m_clearcoat;
    float m_clearcoatGloss;
    float m_specular;
    float m_specTint;
    float m_specTrans;
    float m_roughness;
    float m_anisotropic;
    float m_metallic;
    float m_eta;

    // Secondary parameters
    float m_aspect;
    float m_alpha_min;
    float m_alphax;
    float m_alphay;
    float m_alpha_g;
};

NORI_REGISTER_CLASS(Disney, "disney");
NORI_NAMESPACE_END