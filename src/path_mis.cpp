#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class PathMis : public Integrator {
public:
    PathMis(const PropertyList &props) {

    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const{
        Intersection its;
        Ray3f currRay = ray;
        Color3f t = 1.f;
        Color3f color = 0;
        float w_mat = 1.f;

        while (true) {
            /* return black when no intersection */
            if (!scene->rayIntersect(currRay, its))
                break;

            /* compute Le when intersection is an emitter */
            if (its.mesh->isEmitter()) {
                EmitterQueryRecord lRec(currRay.o, its.p, its.shFrame.n);
                Color3f Li = its.mesh->getEmitter()->eval(lRec);
                color += w_mat * t * its.mesh->getEmitter()->eval(lRec);
            }

            /* sample emitter */
            auto light = scene->getRandomEmitter(sampler->next1D());
            int n_lights = scene->getLights().size();
            EmitterQueryRecord lRec(its.p);
            Color3f Li = light->sample(lRec, sampler->next2D());

            /* compute w_em */
            if (!scene->rayIntersect(lRec.shadowRay)) {  // if shadow ray is NOT occluded
                float cosTheta = its.shFrame.n.dot(lRec.wi);
                BSDFQueryRecord bRec(its.toLocal(-currRay.d), its.toLocal(lRec.wi), ESolidAngle);
                bRec.uv = its.uv;
                float pdf_em = light->pdf(lRec);
                float pdf_mat = its.mesh->getBSDF()->pdf(bRec);
                if (pdf_em + pdf_mat > 1e-8) {
                    float w_em = pdf_em / (pdf_mat + pdf_em);
                    color += w_em * t * its.mesh->getBSDF()->eval(bRec) * Li * cosTheta * n_lights;
                }
            }

            /* Russian roulette */
            float p = std::min(t.maxCoeff(), 0.99f);
            if (sampler->next1D() > p || p <= 0.f) {
                break;
            }
            t /= p;

            /* sample brdf */
            BSDFQueryRecord bRec(its.toLocal((-currRay.d).normalized()));
            bRec.p = its.p;
            bRec.uv = its.uv;
            Color3f brdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());

            /* recursion */
            t *= brdf;
            currRay = Ray3f(bRec.p, its.toWorld(bRec.wo));

            /* current pdf mat */
            float pdf_mat = its.mesh->getBSDF()->pdf(bRec);

            if (!scene->rayIntersect(currRay, its))
                break;

            /* compute w_mat */
            if (its.mesh->isEmitter()) {
                EmitterQueryRecord lRec(currRay.o, its.p, its.shFrame.n);
                float pdf_em = its.mesh->getEmitter()->pdf(lRec);
                if (pdf_em + pdf_mat > 1e-8) {
                    w_mat = pdf_mat / (pdf_mat + pdf_em);
                } else {
                    w_mat = 0.f;
                }
            }

            if (bRec.measure == EDiscrete)
                w_mat = 1.f;
        }

        return color;
    }

    std::string toString() const {
        return "PathMis[]";
    }
};

NORI_REGISTER_CLASS(PathMis, "path_mis")
NORI_NAMESPACE_END