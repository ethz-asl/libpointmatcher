#include "../utest.h"
#include "pointmatcher/OutlierFiltersImpl.h"

using namespace std;
using namespace PointMatcherSupport;

//---------------------------
// Outlier modules
//---------------------------

// Utility classes
class OutlierFilterTest: public IcpHelper
{
public:
	std::shared_ptr<PM::OutlierFilter> testedOutlierFilter;

	// Will be called for every tests
	virtual void SetUp()
	{
		icp.setDefault();
		// Uncomment for console outputs
		//setLogger(PM::get().LoggerRegistrar.create("FileLogger"));
		
		icp.outlierFilters.clear();
	}

	// Will be called for every tests
	virtual void TearDown(){}

	void addFilter(string name, PM::Parameters params)
	{
		testedOutlierFilter = 
			PM::get().OutlierFilterRegistrar.create(name, params);
		icp.outlierFilters.push_back(testedOutlierFilter);
	}

};

//No commun parameters were found for 2D and 3D, tests are spliced
TEST_F(OutlierFilterTest, MaxDistOutlierFilter2D)
{
	addFilter("MaxDistOutlierFilter", {
			{"maxDist", toParam(0.10)}
		}
	);
	validate2dTransformation();
}

TEST_F(OutlierFilterTest, MaxDistOutlierFilter3D)
{
	addFilter("MaxDistOutlierFilter", {
			{"maxDist", toParam(1.0)}
		}
	);
	validate3dTransformation();
}

//No commun parameters were found for 2D and 3D, tests are spliced
TEST_F(OutlierFilterTest, MinDistOutlierFilter2D)
{
	// Since not sure how useful is that filter, we keep the 
	// MaxDistOutlierFilter with it
	std::shared_ptr<PM::OutlierFilter> extraOutlierFilter;
	
	extraOutlierFilter = 
		PM::get().OutlierFilterRegistrar.create(
			"MaxDistOutlierFilter", {
				{"maxDist", toParam(0.10)}
			}
		);
	icp.outlierFilters.push_back(extraOutlierFilter);
	
	addFilter("MinDistOutlierFilter", {{"minDist", toParam(0.0002)}});
	
	validate2dTransformation();
}

TEST_F(OutlierFilterTest, MinDistOutlierFilter3D)
{
	// Since not sure how useful is that filter, we keep the 
	// MaxDistOutlierFilter with it
	std::shared_ptr<PM::OutlierFilter> extraOutlierFilter;
	
	extraOutlierFilter = 
		PM::get().OutlierFilterRegistrar.create(
			"MaxDistOutlierFilter", {
				{"maxDist", toParam(1.0)}
			}
		)
	;
	icp.outlierFilters.push_back(extraOutlierFilter);
	
	addFilter("MinDistOutlierFilter", {{"minDist", toParam(0.0002)}});
	
	validate3dTransformation();
}

TEST_F(OutlierFilterTest, MedianDistOutlierFilter)
{
	addFilter("MedianDistOutlierFilter", {{"factor", toParam(3.5)}});
	validate2dTransformation();
	validate3dTransformation();
}


TEST_F(OutlierFilterTest, TrimmedDistOutlierFilter)
{
	addFilter("TrimmedDistOutlierFilter", {{"ratio", toParam(0.85)}});
	validate2dTransformation();
	validate3dTransformation();
}


TEST_F(OutlierFilterTest, VarTrimmedDistOutlierFilter)
{
	addFilter("VarTrimmedDistOutlierFilter", {
			{"minRatio", toParam(0.60)},
			{"maxRatio", toParam(0.80)},
			{"lambda", toParam(0.9)}
		}
	);
	validate2dTransformation();
	validate3dTransformation();
}

OutlierFiltersImpl<NumericType>::OutlierWeights VarTrimLambdaTest(const NumericType lambda) {
	OutlierFiltersImpl<NumericType>::VarTrimmedDistOutlierFilter filter({{"minRatio", toParam(0.0000001)},
																   {"maxRatio", toParam(1.0)},
																   {"lambda", toParam(lambda)}});
	PointMatcher<NumericType>::DataPoints filteredReading;
	PointMatcher<NumericType>::DataPoints filteredReference;

	// Create a vector a distance
	PointMatcher<NumericType>::Matches::Dists dists(1, 5);
	dists << 4, 5, 5, 5, 5;
	PointMatcher<NumericType>::Matches::Ids ids(1, 5);
	PointMatcher<NumericType>::Matches input(dists, ids);
	return filter.compute(filteredReading, filteredReference, input);
}

TEST_F(OutlierFilterTest, VarTrimmedDistOutlierFilterParameters)
{
	// A lambda parameter of zero, all matches will be reject except for the minimum
	OutlierFiltersImpl<NumericType>::OutlierWeights weights = VarTrimLambdaTest(0.0);
	// The minimum is the first value
	EXPECT_EQ(1.0f, weights(0, 0));
	EXPECT_EQ(0.0f, weights(0, 1));

	weights = VarTrimLambdaTest(1.0);
	EXPECT_EQ(1.0f, weights(0, 0));
	EXPECT_EQ(1.0f, weights(0, 1));
}

// Test for DescriptorMatchOutlierFilter
TEST_F(OutlierFilterTest, DescriptorMatchOutlierFilter)
{
	// Create simple reading and reference point clouds
	PM::Matrix readingFeat(3, 2);
	readingFeat << 0, 1,
		0, 0,
		0, 0;
	PM::Matrix referenceFeat(3, 2);
	referenceFeat << 0, 1.1,
		0, 0,
		0, 0;

	// Define standard feature labels
	PM::DataPoints::Labels featureLabels;
	featureLabels.push_back(PM::DataPoints::Label("x", 1));
	featureLabels.push_back(PM::DataPoints::Label("y", 1));
	featureLabels.push_back(PM::DataPoints::Label("z", 1));

	// Add 'intensity' descriptor
	PM::Matrix readingDesc(1, 2);
	readingDesc << 10, 20; // Intensity values for reading points
	PM::Matrix referenceDesc(1, 2);
	referenceDesc << 10.5, 15; // Intensity values for reference points

	// Create DataPoints with feature labels
	PM::DataPoints reading(readingFeat, featureLabels);
	reading.addDescriptor("intensity", readingDesc);
	PM::DataPoints reference(referenceFeat, featureLabels);
	reference.addDescriptor("intensity", referenceDesc);

	// Create matches (point 0 -> point 0, point 1 -> point 1)
	PM::Matches::Ids ids(1, 2);
	ids << 0, 1;
	PM::Matches::Dists dists(1, 2);
	dists << 0.01, 0.1; // Dummy distances
	PM::Matches matches(dists, ids);

	// --- Test Case: Default (Squared L2 norm, soft threshold) ---
	addFilter("DescriptorMatchOutlierFilter", {
		{"descName", "intensity"},
		{"sigma", toParam(5.0)} // Sigma for weighting
	});

	// Compute weights
	PM::OutlierWeights weights = testedOutlierFilter->compute(reading, reference, matches);

	// Expected weights calculation (Squared L2 soft):
	// Match 0: Reading intensity 10, Reference intensity 10.5. Diff = 0.5. Diff^2 = 0.25
	//          Weight = exp(-0.25 / (5.0^2)) = exp(-0.25 / 25) = exp(-0.01) approx 0.9900
	// Match 1: Reading intensity 20, Reference intensity 15. Diff = 5. Diff^2 = 25
	//          Weight = exp(-25 / (5.0^2)) = exp(-25 / 25) = exp(-1) approx 0.3679

	ASSERT_EQ(weights.rows(), 1);
	ASSERT_EQ(weights.cols(), 2);
	EXPECT_NEAR(weights(0, 0), std::exp(-0.25 / 25.0), 1e-4);
	EXPECT_NEAR(weights(0, 1), std::exp(-25.0 / 25.0), 1e-4);

	// --- Test case with missing descriptor in reference ---
	icp.outlierFilters.clear(); // Use default filter settings for this test
	addFilter("DescriptorMatchOutlierFilter", {{"descName", "intensity"}, {"sigma", toParam(5.0)}});
	PM::DataPoints referenceNoDesc(referenceFeat, featureLabels);
	PM::OutlierWeights weightsMissingRef = testedOutlierFilter->compute(reading, referenceNoDesc, matches);
	EXPECT_EQ(weightsMissingRef(0, 0), 1.0); // Should default to 1 if descriptor missing
	EXPECT_EQ(weightsMissingRef(0, 1), 1.0);
	EXPECT_TRUE(std::dynamic_pointer_cast<OutlierFiltersImpl<NumericType>::DescriptorMatchOutlierFilter>(testedOutlierFilter)->warningPrinted); // Check flag

	// --- Test case with missing descriptor in reading ---
	// Reset warning flag for next test
	std::dynamic_pointer_cast<OutlierFiltersImpl<NumericType>::DescriptorMatchOutlierFilter>(testedOutlierFilter)->warningPrinted = false;
	PM::DataPoints readingNoDesc(readingFeat, featureLabels);
	PM::OutlierWeights weightsMissingRead = testedOutlierFilter->compute(readingNoDesc, reference, matches);
	EXPECT_EQ(weightsMissingRead(0, 0), 1.0); // Should default to 1
	EXPECT_EQ(weightsMissingRead(0, 1), 1.0);
	EXPECT_TRUE(std::dynamic_pointer_cast<OutlierFiltersImpl<NumericType>::DescriptorMatchOutlierFilter>(testedOutlierFilter)->warningPrinted); // Check flag

	// --- Test case with mismatched descriptor dimensions ---
	std::dynamic_pointer_cast<OutlierFiltersImpl<NumericType>::DescriptorMatchOutlierFilter>(testedOutlierFilter)->warningPrinted = false;
	PM::Matrix referenceDesc2D(2, 2);
	referenceDesc2D << 10.5, 15, 1, 2;
	PM::DataPoints referenceMismatchDim(referenceFeat, featureLabels);
	referenceMismatchDim.addDescriptor("intensity", referenceDesc2D);
	PM::OutlierWeights weightsMismatch = testedOutlierFilter->compute(reading, referenceMismatchDim, matches);
	EXPECT_EQ(weightsMismatch(0, 0), 1.0); // Should default to 1
	EXPECT_EQ(weightsMismatch(0, 1), 1.0);
	EXPECT_TRUE(std::dynamic_pointer_cast<OutlierFiltersImpl<NumericType>::DescriptorMatchOutlierFilter>(testedOutlierFilter)->warningPrinted); // Check flag

	// --- Test case with 3D descriptor (RGB) ---
	icp.outlierFilters.clear(); // Clear previous filters
	PM::Matrix readingRgbDesc(3, 2);
	readingRgbDesc << 100, 200, // R
					  50, 150,  // G
					  25, 75;   // B
	PM::Matrix referenceRgbDesc(3, 2);
	referenceRgbDesc << 110, 180, // R
						55, 160,  // G
						30, 80;   // B

	PM::DataPoints readingRgb(readingFeat, featureLabels);
	readingRgb.addDescriptor("rgb", readingRgbDesc);
	PM::DataPoints referenceRgb(referenceFeat, featureLabels);
	referenceRgb.addDescriptor("rgb", referenceRgbDesc);

	addFilter("DescriptorMatchOutlierFilter", {
		{"descName", "rgb"},
		{"sigma", toParam(30.0)} // Sigma for RGB weighting
	});

	PM::OutlierWeights weightsRgb = testedOutlierFilter->compute(readingRgb, referenceRgb, matches);

	// Expected weights calculation (RGB):
	// Match 0: Read [100, 50, 25], Ref [110, 55, 30]. Diff = [-10, -5, -5].
	//          Diff^2 Norm = (-10)^2 + (-5)^2 + (-5)^2 = 100 + 25 + 25 = 150
	//          Weight = exp(-150 / (30.0^2)) = exp(-150 / 900) = exp(-1/6) approx 0.8465
	// Match 1: Read [200, 150, 75], Ref [180, 160, 80]. Diff = [20, -10, -5].
	//          Diff^2 Norm = (20)^2 + (-10)^2 + (-5)^2 = 400 + 100 + 25 = 525
	//          Weight = exp(-525 / (30.0^2)) = exp(-525 / 900) = exp(-7/12) approx 0.5581

	ASSERT_EQ(weightsRgb.rows(), 1);
	ASSERT_EQ(weightsRgb.cols(), 2);
	EXPECT_NEAR(weightsRgb(0, 0), std::exp(-150.0 / 900.0), 1e-4);
	EXPECT_NEAR(weightsRgb(0, 1), std::exp(-525.0 / 900.0), 1e-4);
}

TEST_F(OutlierFilterTest, VarTrimmedDistOutlierFilterWithInfMatches)
{
	OutlierFiltersImpl<NumericType>::VarTrimmedDistOutlierFilter filter({
		{"minRatio", toParam(0.0000001)},
		{"maxRatio", toParam(1.0)},
		{"lambda", toParam(0.0)}
	});
	PointMatcher<NumericType>::DataPoints filteredReading;
	PointMatcher<NumericType>::DataPoints filteredReference;

	// 5 matches: 3 valid finite distances and 2 infinite distances
	const NumericType inf = std::numeric_limits<NumericType>::infinity();
	PointMatcher<NumericType>::Matches::Dists dists(1, 5);
	dists << 2.0, 4.0, inf, 3.0, inf;
	PointMatcher<NumericType>::Matches::Ids ids(1, 5);
	ids << 0, 1, 2, 3, 4;
	PointMatcher<NumericType>::Matches input(dists, ids);

	OutlierFiltersImpl<NumericType>::OutlierWeights weights = filter.compute(filteredReading, filteredReference, input);

	// With lambda=0, only the minimum finite match (distance 2.0 at index 0) should be an inlier
	EXPECT_EQ(1.0, weights(0, 0));
	EXPECT_EQ(0.0, weights(0, 1));
	EXPECT_EQ(0.0, weights(0, 2));
	EXPECT_EQ(0.0, weights(0, 3));
	EXPECT_EQ(0.0, weights(0, 4));
}

TEST_F(OutlierFilterTest, VarTrimmedDistOutlierFilterDegenerateRange)
{
	// 2 matches with minRatio 0.4 and maxRatio 0.45 where floor(0.4*2) == floor(0.45*2) == 0
	OutlierFiltersImpl<NumericType>::VarTrimmedDistOutlierFilter filter({
		{"minRatio", toParam(0.4)},
		{"maxRatio", toParam(0.45)},
		{"lambda", toParam(1.0)}
	});
	PointMatcher<NumericType>::DataPoints filteredReading;
	PointMatcher<NumericType>::DataPoints filteredReference;

	PointMatcher<NumericType>::Matches::Dists dists(1, 2);
	dists << 1.0, 5.0;
	PointMatcher<NumericType>::Matches::Ids ids(1, 2);
	ids << 0, 1;
	PointMatcher<NumericType>::Matches input(dists, ids);

	EXPECT_NO_THROW({
		OutlierFiltersImpl<NumericType>::OutlierWeights weights = filter.compute(filteredReading, filteredReference, input);
		EXPECT_EQ(1.0, weights(0, 0));
	});
}

TEST_F(OutlierFilterTest, VarTrimmedDistOutlierFilterFloatPrecision)
{
	OutlierFiltersImpl<float>::VarTrimmedDistOutlierFilter filter({
		{"minRatio", toParam(0.1)},
		{"maxRatio", toParam(0.9)},
		{"lambda", toParam(1.5)}
	});
	PointMatcher<float>::DataPoints filteredReading;
	PointMatcher<float>::DataPoints filteredReference;

	PointMatcher<float>::Matches::Dists dists(1, 10);
	dists << 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f;
	PointMatcher<float>::Matches::Ids ids(1, 10);
	ids.setZero();
	PointMatcher<float>::Matches input(dists, ids);

	EXPECT_NO_THROW({
		OutlierFiltersImpl<float>::OutlierWeights weights = filter.compute(filteredReading, filteredReference, input);
		EXPECT_EQ(weights.rows(), 1);
		EXPECT_EQ(weights.cols(), 10);
	});
}
