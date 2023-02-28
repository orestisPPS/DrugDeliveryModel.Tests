using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.DrugDeliveryModel.Tests.EquationModels;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.Staggered;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.Constitutive.Structural;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.NumericalAnalyzers;
using MGroup.Solvers.Direct;
using Xunit;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using MGroup.Constitutive.ConvectionDiffusion.InitialConditions;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Solution;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class EquationTCell
    {
        //----------------------------Equation T - Tumor Cell Density-------------------------------
        
        //------------------------------------Variables---------------------------------------------
        /// <summary>
        /// The average value of the three components of the solid tumor velocity vector  [m/s]
        /// </summary>
        private const double SolidSpeed = 2.32E-4; // [m/s]

        /// <summary>
        /// Growth rate parameter 1[mol/(m3)]
        /// </summary>
        private const double K1 = 1.74E-6; // [1/s]
        
        /// <summary>
        /// Growth rate parameter 2[mol/(m3)]
        /// </summary>
        private const double K2 = 8.3E-3; // [mol/(m3)]
        
        /// <summary>
        /// Oxygen concentration (Dependent Variable) [mol/m3]
        /// </summary>
        private const double Cox = 0.2; // [mol / m3]
            
        //---------------------------------------Logging----------------------------------
        /// <summary>
        /// The degree of freedom that will be monitored for equation cox
        /// </summary>
        private ConvectionDiffusionDof coxMonitorDOF = ConvectionDiffusionDof.UnknownVariable;

        /// <summary>
        /// The coordinates of the monitored node
        /// </summary>
        //private double[] monitorNodeCoords = { 0.05, 0.05, 0.05 };
        private double[] monitorNodeCoords = { 0.0, 0.0, 0.0 };
        

        //---------------------------------------Time Discretization Specs------------------------------
        private const double TotalTime = 1E-2;
        
        /// <summary>
        /// For increased accuracy use time-step of order 1E-5
        /// </summary>
        private const double TimeStep = 1E-5;// sec
        
        /// <summary>
        /// Simplified version of the production term without non-linear term
        /// </summary>
        //private readonly Func<double> dependantSourceCoefficient =() => (K1 * Cox) / (K2 + Cox);
        private readonly Func<double> dependantSourceCoefficient =() => 0d;
        
        //---------------------------------------Initial conditions------------------------------
        private double initialTCellDensity = 500d;
        
        private Model model;

        private ComsolMeshReader mesh;

        private List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> neumannBC
            = new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>();

        private static ConvectionDiffusionDof[] constrainedDofType = new ConvectionDiffusionDof[1] { ConvectionDiffusionDof.UnknownVariable };
        private static double[] boundaryValue = new double[1] { 0d };
        private static List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> dirichletBC =
            new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>()
            {
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.TopDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0.1}}, boundaryValue),
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.RightDirichlet, constrainedDofType, new double[1][]{new double[3] {0.1,0,0}}, boundaryValue),
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.FrontDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0}}, boundaryValue),
            };

        public EquationTCell()
        {
            
        }
        
        [Theory]
        [InlineData("../../../DataFiles/workingTetMesh2185_1Domain.mphtxt")]
        public void SolveEquationTc(string fileName)
        {
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
            
            mesh = new ComsolMeshReader(fileName);
            
            model = GetModel();
            
            BoundaryAndInitialConditionsUtility.AssignConvectionDiffusionDirichletBCsToModel(model, dirichletBC, 1E-5);
            BoundaryAndInitialConditionsUtility.AssignConvectionDiffusionICToModel(model, initialTCellDensity);
            
            var dynamicAnalyzer = GetAppropriateSolverAnalyzerAndLog(model).Item1;
            dynamicAnalyzer.Solve();

            var totalNewmarkStepsNum = (int)Math.Truncate(TotalTime / TimeStep);
            var tCell = new double[totalNewmarkStepsNum];
            for (int i1 = 0; i1 < totalNewmarkStepsNum; i1++)
            {
                var timeStepResultsLog = dynamicAnalyzer.Logs[i1];
                tCell[i1] = ((DOFSLog)timeStepResultsLog).DOFValues[model.GetNode(GetAppropriateSolverAnalyzerAndLog(model).Item2), coxMonitorDOF];
            }
            
            CSVExporter.ExportVectorToCSV(tCell, "../../../Integration/Tc_nodes_mslv.csv");
        }

        
        
        private Model GetModel()
        {
            var capacity = 1;
            var diffusionCoefficient = 0d;
            var convectionCoefficient = SolidSpeed;
            var dependentProductionCoefficient = dependantSourceCoefficient();
            var independentSourceCoefficient = 0d;
            
            //Assign equation properties to the domain elements
            var convectionDomainCoefficients = new Dictionary<int, double[]>();   
            var dependentProductionCoefficients = new Dictionary<int, double>();
            var independentProductionCoefficients = new Dictionary<int, double>();
            
            foreach (var elementConnectivity in mesh.ElementConnectivity)
            {
                convectionDomainCoefficients[elementConnectivity.Key] = new double[] { convectionCoefficient, convectionCoefficient, convectionCoefficient };
                dependentProductionCoefficients[elementConnectivity.Key] = dependentProductionCoefficient;
                independentProductionCoefficients[elementConnectivity.Key] = independentSourceCoefficient;
            }

            //Create Model
            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(mesh);
            var model = modelProvider.CreateModelFromComsolFile(convectionDomainCoefficients, diffusionCoefficient, 
                dependentProductionCoefficients, independentProductionCoefficients, capacity);
            
            return model;
        }
        
        
        private (IParentAnalyzer, int)  GetAppropriateSolverAnalyzerAndLog(Model model)
        {
            
            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false };
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var problem = new ProblemConvectionDiffusion(model, algebraicModel);

            var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);

            var dynamicAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, problem, linearAnalyzer, timeStep: TimeStep, totalTime: TotalTime); 
            //var dynamicAnalyzerBuilder = new BDFDynamicAnalyzer.Builder(algebraicModel, problem, linearAnalyzer, timeStep: TimeStep, totalTime: TotalTime, bdfOrder: 5);
            var dynamicAnalyzer = dynamicAnalyzerBuilder.Build();

            // Create a log for the desired dof
            var nodeIdToMonitor =Utilities.FindNodeIdFromNodalCoordinates(mesh.NodesDictionary, monitorNodeCoords, 1e-3);
            var watchDofs = new List<(INode node, IDofType dof)>()
            {
                (model.NodesDictionary[nodeIdToMonitor], coxMonitorDOF),
            };

            linearAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs, algebraicModel);
            dynamicAnalyzer.ResultStorage = new ImplicitIntegrationAnalyzerLog();
            dynamicAnalyzer.Initialize();
            return (dynamicAnalyzer, nodeIdToMonitor);
        }
        
    }
}
