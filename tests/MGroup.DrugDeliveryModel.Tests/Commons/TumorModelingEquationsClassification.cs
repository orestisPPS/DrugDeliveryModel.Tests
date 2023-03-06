namespace MGroup.DrugDeliveryModel.Tests.Commons;

/// <summary>
/// Tumor Growth Model Equations.
/// </summary>
enum TumorGrowthModelEquations
{
    StructuralU,
    FluidPressurePi,
    OxygenConcentrationCOx,
    CancerCellDensityT,
    GrowthRatioLg,
}

public enum TumorGrowthTestCases
{
    /// <summary>
    /// One way coupling: Solid phase contributes to fluid phase.
    /// Pressure Dirichlet : pi=0 -> Left, Right, Front, Back, Top, Bottom.
    /// Structural Dirichlet : u * n = 0 (roller) -> Left,Front, Back, Bottom.
    /// Structural Neumann : 4 Point Loads at the 4 corner nodes of the Top Face.
    /// </summary>
    StructuralToPressureTopRightFree4PointLoadsTop,
    
    /// <summary>
    /// One way coupling: Solid phase contributes to fluid phase.
    /// Pressure Dirichlet : pi=0 -> Left, Right, Front, Back, Top, Bottom.
    /// Structural Dirichlet : u * n = 0 (roller) -> Left, Right, Front, Back, Bottom.
    /// Structural Neumann : 4 Point Loads at the 4 corner nodes of the Top Face.
    /// </summary>
    StructuralToPressureTopFree4PointLoadsTop,
    
    /// <summary>
    /// One way coupling: Fluid phase contributes to solid phase.
    /// Simplified Source term.
    /// Pressure Dirichlet : pi=0 -> Left, Right, Front, Back, Top, Bottom.
    /// Structural Dirichlet : u * n = 0 (roller) -> Left,Front, Back, Bottom.
    /// Structural Neumann : Zero Flux -> Top
    /// </summary>
    PressureToStructural,
    
    /// <summary>
    /// Full coupling: Fluid phase contributes to solid phase and vice versa.
    /// Simplified Source term.
    /// Pressure Dirichlet : pi=0 -> Left, Right, Front, Back, Top, Bottom.
    /// Structural Dirichlet : u * n = 0 (roller) -> Left,Front, Back, Bottom.
    /// Structural Neumann : 4 Point Loads at the 4 corner nodes of the Top Face.
    /// </summary>
    StructuralPressureCoupledSimplifiedSourceTerm,
    
    /// <summary>
    /// Full coupling: Fluid phase contributes to solid phase and vice versa.
    /// Complete Source term.
    /// Pressure Dirichlet : pi=0 -> Left, Right, Front, Back, Top, Bottom.
    /// Structural Dirichlet : u * n = 0 (roller) -> Left,Front, Back, Bottom.
    /// Structural Neumann : 4 Point Loads at the 4 corner nodes of the Top Face.
    /// </summary>
    StructuralPressureCoupled,
    
}

	public enum BoundaryConditionCase
	{
		LeftDirichlet,
		LeftPointFlux,
		LeftDistributedFlux,
		
		RightDirichlet,
		RightPointFlux,
		RightDistributedFlux,
		
		FrontDirichlet,
		FrontPointFlux,
		FrontDistributedFlux,
		
		BackDirichlet,
		BackPointFlux,
		BackDistributedFlux,
		
		BottomDirichlet,
		BottomPointFlux,
		BottomDistributedFlux,
		
		TopDirichlet,
		TopPointFlux,
		TopDistributedFlux
	}
	
	

public class TumorModelingEquationsClassification
{

}