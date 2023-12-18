#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/phasefunction.h>
#include <stack>
NORI_NAMESPACE_BEGIN
using namespace std;

class VolPathMATSIntegrator : public Integrator {
private:
    float ray_length;
public:
    VolPathMATSIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Color3f color = 0;
        // larger throughput => larger contribution => larger probability for russian
        Color3f t = 1.f;
        Ray3f currRay = ray;
        int bounce_cnt = 0;
        const Medium* in_medium = nullptr;
        Intersection its;
        stack<Medium> media;
        bool has_intersection = scene->rayIntersect(currRay, its);
        while (true) {
            MediumQueryRecord mRec(currRay.o, -currRay.d, its.t);
            if (in_medium && in_medium->sample_freepath(mRec, sampler)) {
                // scattered inside the medium
                PhaseFunctionQueryRecord pRec(mRec.wi);
                in_medium->getPhaseFunction()->sample(pRec, sampler->next2D());   //sample direction to next interaction
                
                currRay = Ray3f(mRec.p, pRec.wo);
                has_intersection = scene->rayIntersect(currRay, its);
                t *= mRec.albedo;
            }
            else {
                if (in_medium) t *= mRec.albedo;
                if (!has_intersection) {
                    break;
                }
                if (its.mesh->isEmitter()) {
                    EmitterQueryRecord lRec(currRay.o, its.p, its.shFrame.n);
                    color += t * its.mesh->getEmitter()->eval(lRec);
                }
                // sample BSDF ray

                BSDFQueryRecord bRec(its.shFrame.toLocal(-currRay.d));
                bRec.uv = its.uv;
                bRec.p = its.p;
                Color3f brdf= its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                currRay = Ray3f(bRec.p, its.shFrame.toWorld(bRec.wo));
                in_medium = its.mesh->getMedium();  
                has_intersection = scene->rayIntersect(currRay, its);
                //in_medium = its.mesh->getMedium();
                t *= brdf;
            }
            // Ruassian Roulette 
            if (++bounce_cnt > 3) {
                float p = std::min(t.maxCoeff(), 0.99f);
                if (sampler->next1D() > p ) break;
                t /= p;
            }
        }
        return color;
    }

    std::string toString() const {
        return "VolPathMATSIntegrator[]";
    }
};

NORI_REGISTER_CLASS(VolPathMATSIntegrator, "vol_path_mats");
NORI_NAMESPACE_END