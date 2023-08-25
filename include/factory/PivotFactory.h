//
// Created by joaovictor on 17/08/23.
//

#ifndef GERVLIB_PIVOTFACTORY_H
#define GERVLIB_PIVOTFACTORY_H

#include "BPPPivots.h"
#include "ConvexPivots.h"
#include "FFTPivots.h"
#include "HFIPivots.h"
#include "ISPivots.h"
#include "KmedoidsPivots.h"
#include "MaxSeparatedPivots.h"
#include "MaxVariancePivots.h"
#include "PCAPivots.h"
#include "RandomPivots.h"
#include "SelectionPivots.h"
#include "SSSPivots.h"
#include "WDRPivots.h"

namespace gervLib::pivots
{

    template <typename O, typename T>
    class PivotFactory
    {

    public:
        static std::unique_ptr<pivots::Pivot<O, T>> createPivot(pivots::PIVOT_TYPE type)
        {

            if (type == RANDOM)
                return std::make_unique<RandomPivots<O, T>>();
            else if (type == SELECTION)
                return std::make_unique<SelectionPivots<O, T>>();
            else if (type == KMEDOIDS)
                return std::make_unique<KmedoidsPivots<O, T>>();
            else if (type == PCA)
                return std::make_unique<PCAPivots<O, T>>();
            else if (type == FFT)
                return std::make_unique<FFTPivots<O, T>>();
            else if (type == BPP)
                return std::make_unique<BPPPivots<O, T>>();
            else if (type == HFI)
                return std::make_unique<HFIPivots<O, T>>();
            else if (type == WDR)
                return std::make_unique<WDRPivots<O, T>>();
            else if (type == MAXVARIANCE)
                return std::make_unique<MaxVariancePivots<O, T>>();
            else if (type == MAXSEPARATED)
                return std::make_unique<MaxSeparatedPivots<O, T>>();
            else if (type == CONVEX)
                return std::make_unique<ConvexPivots<O, T>>();
            else if (type == IS)
                return std::make_unique<ISPivots<O, T>>();
            else if (type == SSS)
                return std::make_unique<SSSPivots<O, T>>();
            else
                throw std::invalid_argument("Invalid PIVOT_TYPE");
        }

        static std::unique_ptr<pivots::Pivot<O, T>> createPivot(pivots::PIVOT_TYPE type, std::unique_ptr<pivots::Pivot<O, T>>& pvt)
        {

            std::unique_ptr<pivots::Pivot<O, T>> ans;

            if (type == RANDOM)
                ans = std::make_unique<RandomPivots<O, T>>();
            else if (type == SELECTION)
                ans = std::make_unique<SelectionPivots<O, T>>();
            else if (type == KMEDOIDS)
                ans = std::make_unique<KmedoidsPivots<O, T>>();
            else if (type == PCA)
                ans = std::make_unique<PCAPivots<O, T>>();
            else if (type == FFT)
                ans = std::make_unique<FFTPivots<O, T>>();
            else if (type == BPP)
                ans = std::make_unique<BPPPivots<O, T>>();
            else if (type == HFI)
                ans = std::make_unique<HFIPivots<O, T>>();
            else if (type == WDR)
                ans = std::make_unique<WDRPivots<O, T>>();
            else if (type == MAXVARIANCE)
                ans = std::make_unique<MaxVariancePivots<O, T>>();
            else if (type == MAXSEPARATED)
                ans = std::make_unique<MaxSeparatedPivots<O, T>>();
            else if (type == CONVEX)
                ans = std::make_unique<ConvexPivots<O, T>>();
            else if (type == IS)
                ans = std::make_unique<ISPivots<O, T>>();
            else if (type == SSS)
                ans = std::make_unique<SSSPivots<O, T>>();
            else
                throw std::invalid_argument("Invalid PIVOT_TYPE");

            ans->setSeed(pvt->getSeed());
            ans->setSampleSize(pvt->getSampleSize());
            ans->setNumberOfDropPivots(pvt->getNumberOfDropPivots());

            return ans;

        }

    };

}


#endif //GERVLIB_PIVOTFACTORY_H
