#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class PathMats : public Integrator {
public:
    PathMats(const PropertyList &props) {

    }
        
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;
        Ray3f currRay = ray;
        Color3f t = 1.f;
        Color3f color = 0;

        /* path tracing */
        while (true) {
            if (!scene->rayIntersect(currRay, its)) {
                break;
            }

            /* compute Le when intersection is an emitter */
            if (its.mesh->isEmitter()) {
                EmitterQueryRecord lRec(currRay.o, its.p, its.shFrame.n);
                color += t * its.mesh->getEmitter()->eval(lRec);
            }

            /* Russian roulette */
            float p = std::min(t.maxCoeff(), 0.99f);
            if (sampler->next1D() > p || p <= 0.f) {
                break;
            }
            t /= p;

            /* sample brdf */
            BSDFQueryRecord bRec(its.shFrame.toLocal(-currRay.d));
            bRec.p = its.p;
            bRec.uv = its.uv;
            Color3f brdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            t *= brdf;

            /* recursion */
            currRay = Ray3f(bRec.p, its.toWorld(bRec.wo));
        }

        return color;
    }

    std::string toString() const {
        return "PathMats[]";
    }
};

NORI_REGISTER_CLASS(PathMats, "path_mats")
NORI_NAMESPACE_END