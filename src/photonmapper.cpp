/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/scene.h>
#include <nori/photon.h>

NORI_NAMESPACE_BEGIN

class PhotonMapper : public Integrator {
public:
    /// Photon map data structure
    typedef PointKDTree<Photon> PhotonMap;
    int m_emittedPhotonCount = 0;

    PhotonMapper(const PropertyList &props) {
        /* Lookup parameters */
        m_photonCount  = props.getInteger("photonCount", 1000000);
        m_photonRadius = props.getFloat("photonRadius", 0.0f /* Default: automatic */);
    }

    virtual void preprocess(const Scene *scene) override {
        cout << "Gathering " << m_photonCount << " photons .. ";
        cout.flush();

        /* Create a sample generator for the preprocess step */
        Sampler *sampler = static_cast<Sampler *>(
            NoriObjectFactory::createInstance("independent", PropertyList()));

        /* Allocate memory for the photon map */
        m_photonMap = std::unique_ptr<PhotonMap>(new PhotonMap());
        m_photonMap->reserve(m_photonCount);

		/* Estimate a default photon radius */
		if (m_photonRadius == 0)
			m_photonRadius = scene->getBoundingBox().getExtents().norm() / 500.0f;
	

		/* How to add a photon?
		 * m_photonMap->push_back(Photon(
		 *	Point3f(0, 0, 0),  // Position
		 *	Vector3f(0, 0, 1), // Direction
		 *	Color3f(1, 2, 3)   // Power
		 * ));
		 */

		// put your code to trace photons here

        /* trace photon */
        int currPhotonCount = 0;
        m_emittedPhotonCount = 0;
        while (currPhotonCount < m_photonCount) {
            /* sample photon */
            const Emitter* light = scene->getRandomEmitter(sampler->next1D());
            const int n_lights = scene->getLights().size();
            Ray3f currRay;
            Color3f t = 1.f;
            Color3f w = light->samplePhoton(currRay, sampler->next2D(), sampler->next2D()) * n_lights;
            Intersection its;
            m_emittedPhotonCount ++;
            
            /* trace photon */
            while (true) {
                /* return black when no intersection */
                if (!scene->rayIntersect(currRay, its))
                    break;

                /* Russian roulette */
                float p = std::min(t.maxCoeff(), 0.99f);
                if (sampler->next1D() > p || p <= 0.f) {
                    break;
                }
                t /= p;

                /* if diffuse surface, add photon record */
                if (its.mesh->getBSDF()->isDiffuse()) {
                    m_photonMap->push_back(Photon(
                        its.p,
                        -currRay.d,
                        t * w
                    ));
                    currPhotonCount ++;
                    if (currPhotonCount == m_photonCount) {
                        break;
                    }
                }

                /* sample brdf */
                BSDFQueryRecord bRec(its.toLocal(-currRay.d));
                bRec.p = its.p;
                bRec.uv = its.uv;
                Color3f brdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());

                /* recursion */
                t *= brdf;
                currRay = Ray3f(bRec.p, its.toWorld(bRec.wo));
            }
        }

		/* Build the photon map */
        m_photonMap->build();
    }

    virtual Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const override {
    	
		/* How to find photons?
		 * std::vector<uint32_t> results;
		 * m_photonMap->search(Point3f(0, 0, 0), // lookup position
		 *                     m_photonRadius,   // search radius
		 *                     results);
		 *
		 * for (uint32_t i : results) {
		 *    const Photon &photon = (*m_photonMap)[i];
		 *    cout << "Found photon!" << endl;
		 *    cout << " Position  : " << photon.getPosition().toString() << endl;
		 *    cout << " Power     : " << photon.getPower().toString() << endl;
		 *    cout << " Direction : " << photon.getDirection().toString() << endl;
		 * }
		 */

		// put your code for path tracing with photon gathering here
        Intersection its;
        Ray3f currRay = _ray;
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

            /* if diffuse surface, estimate photon density */
            if (its.mesh->getBSDF()->isDiffuse()) {
                Color3f Lp = 0;

                /* find photons */
                std::vector<uint32_t> results;
		        m_photonMap->search(its.p, // lookup position
		            m_photonRadius,   // search radius
		            results);
                
                /* photon gathering */
                float area = M_PI * m_photonRadius * m_photonRadius;
                for (uint32_t i : results) {
                    const Photon &photon = (*m_photonMap)[i];
                    BSDFQueryRecord bRec(its.shFrame.toLocal(-currRay.d), its.shFrame.toLocal(photon.getDirection()), ESolidAngle);
                    bRec.p = its.p;
                    bRec.uv = its.uv;
                    Lp += its.mesh->getBSDF()->eval(bRec) * photon.getPower() / m_emittedPhotonCount;
                }

                color += t * (Lp * INV_PI / (m_photonRadius * m_photonRadius));
                break;
            }

            /* Russian roulette */
            float p = std::min(t.maxCoeff(), 0.99f);
            // float p = 0.95;
            if (sampler->next1D() > p || p <= 0.f) {
                break;
            }
            t /= p;

            /* sample brdf */
            BSDFQueryRecord bRec(its.toLocal(-currRay.d));
            bRec.p = its.p;
            bRec.uv = its.uv;
            Color3f brdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            t *= brdf;

            /* recursion */
            currRay = Ray3f(bRec.p, its.toWorld(bRec.wo));
        }

		return color;
    }

    virtual std::string toString() const override {
        return tfm::format(
            "PhotonMapper[\n"
            "  photonCount = %i,\n"
            "  photonRadius = %f\n"
            "]",
            m_photonCount,
            m_photonRadius
        );
    }
private:
    /* 
     * Important: m_photonCount is the total number of photons deposited in the photon map,
     * NOT the number of emitted photons. You will need to keep track of those yourself.
     */ 
    int m_photonCount;
    float m_photonRadius;
    std::unique_ptr<PhotonMap> m_photonMap;
};

NORI_REGISTER_CLASS(PhotonMapper, "photonmapper");
NORI_NAMESPACE_END
