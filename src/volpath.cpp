#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/phasefunction.h>
#include <stack>
NORI_NAMESPACE_BEGIN
using namespace std;

class VolPathMISIntegrator : public Integrator {
private:
    float ray_length;
public:
    VolPathMISIntegrator(const PropertyList& props) {
    }
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Color3f color = 0;
        Color3f t = 1.f;
        Ray3f currRay = ray;
        Ray3f trackRay;
        int bounce_cnt = 0;
        const Medium* current_medium = nullptr;
        Intersection its;
        Intersection trackIts;
        stack<const Medium*> media;
        float w_mat = 1.f;
        bool has_intersection = scene->rayIntersect(currRay, its);
        current_medium = scene->getCameraMedium();
        while (true) {
            MediumQueryRecord mRec(currRay.o, -currRay.d, its.t);

            if (current_medium && current_medium->sample_freepath(mRec, sampler)) {
                // continuously scattering inside the medium
                PhaseFunctionQueryRecord pRec(mRec.wi);
                current_medium->getPhaseFunction()->sample(pRec, sampler->next2D());   
                //sample direction to next interaction
                currRay = Ray3f(mRec.p, pRec.wo);
                has_intersection = scene->rayIntersect(currRay, its);
                t *= mRec.ret;
            }
            else {
                if (current_medium) {
                    //hit a surface of mesh or medium during scattering, update the throughput
                    t *= mRec.ret;
                }
                if (!has_intersection) {
                    break;
                }

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
                if (Frame::cosTheta(its.shFrame.toLocal(currRay.d)) >= 0.0f && its.mesh->getMedium() == current_medium) {
                    //if escaping the medium, carefully pop
                    if (!media.empty()) {
                        media.pop();
                        if (media.empty()) {
                            current_medium = nullptr;
                        }
                        else {
                            current_medium = media.top();
                        }
                    }
                }
                else if (Frame::cosTheta(its.shFrame.toLocal(currRay.d)) < 0.0f && its.mesh->getMedium() != nullptr) {
                    //if entering medium, push into the stack
                    current_medium = its.mesh->getMedium();
                    media.push(current_medium);
                }

                has_intersection = scene->rayIntersect(currRay, its);

                /* compute w_mat */
                if (its.mesh->isEmitter()) {
                    EmitterQueryRecord lRec(currRay.o, its.p, its.shFrame.n);
                    float pdf_em = its.mesh->getEmitter()->pdf(lRec);
                    if (pdf_em + pdf_mat > 1e-8) {
                        w_mat = pdf_mat / (pdf_mat + pdf_em);
                    }
                    else {
                        w_mat = 0.f;
                    }
                }

                if (bRec.measure == EDiscrete)
                    w_mat = 1.f;
            }
        }

        return color;
    }

    std::string toString() const {
        return "VolPathMISIntegrator[]";
    }
};

NORI_REGISTER_CLASS(VolPathMISIntegrator, "vol_path_mis");
NORI_NAMESPACE_END