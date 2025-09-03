template<class DataPoint, class _WFunctor, typename T>
template <typename IndexRange, typename PointContainer>
PONCA_MULTIARCH inline FIT_RESULT MLS<DataPoint, _WFunctor, T>::computeWithIdsMLS(IndexRange ids, const PointContainer& points){
    lastPosMLS = Base::getWeightFunc().basisCenter(); // Initial position before projection
    FIT_RESULT res = FIT_RESULT::UNDEFINED;

    for( int mm = 0; mm < mlsIter; ++mm) {
        // Starts a new pass and initialise the fit
        Base::computeWithIds(ids, points); // TODO : Fix this unresolved call
        res = Base::finalize();

        if (Base::isStable()) {
            lastPosMLS = Base::project( lastPosMLS );
        } else {
            std::cerr << "Warning: fit at mls iteration " << mm << " is not stable" << std::endl;
            break;
        }
    }
    return res;
};
