# trace generated using paraview version 5.10.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'OpenFOAMReader'
rundx50foam = OpenFOAMReader(registrationName='rundx50.foam', FileName='/mnt/e/CFPremixed-DMEAir/DME-phi_01.00/GridConvergence/DMEIgn_Tu600_U1.0_L20.0S_R30.0_R500_T200_dx50/rundx50.foam')

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on rundx50foam
rundx50foam.Adddimensionalunitstoarraynames = 1

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=rundx50foam)

UpdatePipeline(time=5e-05, proxy=plotOverLine1)

# Properties modified on plotOverLine1
plotOverLine1.Point1 = [0.0, 0.0, 0.0]
plotOverLine1.Point2 = [0.01, 0.0, 0.0]

# Properties modified on plotOverLine1
plotOverLine1.SamplingPattern = 'Sample At Segment Centers'

# create extractor
cSV1 = CreateExtractor('CSV', plotOverLine1, registrationName='CSV1')
# trace defaults for the extractor.
cSV1.Trigger = 'TimeStep'

# Properties modified on cSV1.Writer
cSV1.Writer.FileName = 'lineY0.0_{time:06f}.csv'

# set active source
SetActiveSource(plotOverLine1)

# set active source
SetActiveSource(cSV1)

# save extracts
SaveExtracts(ExtractsOutputDirectory='/mnt/e/CFPremixed-DMEAir/DME-phi_01.00/GridConvergence/DMEIgn_Tu600_U1.0_L20.0S_R30.0_R500_T200_dx50/S3.0E09/ParaViewExtract/lineY0.0',
    GenerateCinemaSpecification=0)

# set active source
SetActiveSource(rundx50foam)

# create a new 'Plot Over Line'
plotOverLine2 = PlotOverLine(registrationName='PlotOverLine2', Input=rundx50foam)

UpdatePipeline(time=0.003, proxy=plotOverLine2)

# Properties modified on plotOverLine2
plotOverLine2.Point1 = [0.0, 0.0, 0.0]
plotOverLine2.Point2 = [0.0, 0.03, 0.0]

# Properties modified on plotOverLine2
plotOverLine2.SamplingPattern = 'Sample At Segment Centers'

# create extractor
cSV2 = CreateExtractor('CSV', plotOverLine2, registrationName='CSV2')

# Properties modified on cSV2
cSV2.Trigger = 'TimeStep'

# Properties modified on cSV2.Writer
cSV2.Writer.FileName = 'lineX0.0_{time:.06f}.csv'

# save extracts
SaveExtracts(ExtractsOutputDirectory='/mnt/r/CFPremixed-DMEAir/DME-phi_01.00/GridConvergence/DMEIgn_Tu600_U1.0_L20.0S_R30.0_R500_T200_dx50/S3.0E09/ParaViewExtract/lineX0.0',
    GenerateCinemaSpecification=0)
