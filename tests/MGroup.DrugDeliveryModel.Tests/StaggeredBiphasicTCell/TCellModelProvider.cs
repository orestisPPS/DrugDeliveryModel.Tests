using System;
using System.Collections.Generic;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.DrugDeliveryModel.Tests.EquationModels;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.Solvers.Direct;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Solution;
using MGroup.NumericalAnalyzers;

namespace MGroup.DrugDeliveryModel.Tests.PreliminaryModels;

public class TCellModelProvider
{
    private double K1 { get; }
    private double K2 { get; }
    
    private Dictionary<int, double> DomainCOx { get; }
    
    private Dictionary<int, double[]> SolidVelocityDivergence { get; }
    private ComsolMeshReader Mesh { get; }
    private ConvectionDiffusionDof MonitorDOFType { get; }
    private int MonitorNodeId { get; }
    private List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> DirichletBCs { get; }
    private List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> NeumannBCs { get; }
    
    private double InitialCondition { get; }

    public TCellModelProvider(double k1, double k2, Dictionary<int, double> domainCOx, Dictionary<int, double[]> solidVelocityDivergence,
        ComsolMeshReader mesh,
        ConvectionDiffusionDof tCellMonitorDOFType, int monitorNodeId,
        List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> dirichletBCs,
        List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> neumannBCs,
        double initialCondition)
    {
        K1 = k1;
        K2 = k2;
        DomainCOx = domainCOx;
        SolidVelocityDivergence = solidVelocityDivergence;
        Mesh = mesh;
        MonitorDOFType = tCellMonitorDOFType;
        MonitorNodeId = monitorNodeId;
        DirichletBCs = dirichletBCs;
        NeumannBCs = neumannBCs;
        InitialCondition = initialCondition;

        IsoparametricJacobian3D.DeterminantTolerance = 1E-20;
    }
    
    public Model GetModel()
    {
        var capacity = 1;
        var diffusionCoefficient = 0d;

            
        //Assign equation properties to the domain elements
        var convectionDomainCoefficients = new Dictionary<int, double[]>();   
        var dependentProductionCoefficients = new Dictionary<int, double>();
        var independentProductionCoefficients = new Dictionary<int, double>();
            
        foreach (var elementConnectivity in Mesh.ElementConnectivity)
        {
            var vs = SolidVelocityDivergence[elementConnectivity.Key];
            convectionDomainCoefficients[elementConnectivity.Key] = new double [] {vs[0], vs[0], vs[0]};
            
            var elementCOx = DomainCOx[elementConnectivity.Key];
            var dependentProductionCoefficient = (K1 * elementCOx) / (K2 + elementCOx);
            dependentProductionCoefficients[elementConnectivity.Key] = dependentProductionCoefficient;
            
            independentProductionCoefficients[elementConnectivity.Key] = 0d;
        }

        //Create Model
        var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(Mesh);
        var model = modelProvider.CreateModelFromComsolFile(convectionDomainCoefficients, diffusionCoefficient, 
            dependentProductionCoefficients, independentProductionCoefficients, capacity);
            
        return model;
    }
    
    public void AddBoundaryConditions(Model model)
    {
        BoundaryAndInitialConditionsUtility.AssignConvectionDiffusionDirichletBCsToModel(model, DirichletBCs, 1E-5);
    }
    
    public void AddInitialConditions(Model model)
    {
        BoundaryAndInitialConditionsUtility.AssignConvectionDiffusionICToModel(model, InitialCondition);
    }
    
        public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog
        (Model model, double pseudoTimeStep, double pseudoTotalTime, int currentStep)
        {
            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false }; //Dense Matrix Solver solves with zero matrices!
            //var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var provider = new ProblemConvectionDiffusion(model, algebraicModel);

    
           var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, provider);

            var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, provider, numIncrements: 1)
            {
                ResidualTolerance = 1E-8,
                MaxIterationsPerIncrement = 100,
                NumIterationsForMatrixRebuild = 1
            };
            var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();

            loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(new List<(INode node, IDofType dof)>()
            {(model.NodesDictionary[MonitorNodeId], MonitorDOFType)}, algebraicModel);

            var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, false, currentStep: currentStep);
            analyzerBuilder.SetNewmarkParametersForConstantAcceleration();
            var analyzer = analyzerBuilder.Build();
            var watchDofs = new[]
            {
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[MonitorNodeId], MonitorDOFType),
                }
            };
            loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);

             return (analyzer, solver, loadControlAnalyzer);
        }
}