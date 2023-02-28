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


        public EquationTCell()
        {
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }
        
        [Theory]
        [InlineData("../../../DataFiles/workingTetMesh2185_1Domain.mphtxt")]
        public void SolveEquationTc(string fileName)
        {
            var capacity = 1;
            var diffusionCoefficient = 0d;
            var convectionCoefficient = SolidSpeed;
            var dependentProductionCoefficient = dependantSourceCoefficient();
            var independentSourceCoefficient = 0d;

            //Read Mesh From comsol file
            var mesh = new ComsolMeshReader(fileName);
            
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

            //Assign Boundary Conditions
            //AddTopRightBackNodesBC(model, 0d, 0, 0.1, 0, 0.1, 0, 0.1);
            //AddInitialConditions(model, 500d, 0, 0.1, 0, 0.1, 0, 0.1);
            
            AddTopRightBackNodesBC(model, 500d, 0, 0.1, 0, 0.1, 0, 0.1);
            AddInitialConditions(model, 50d, 0, 0.1, 0, 0.1, 0, 0.1);
            
            //AddTopRightBackNodesBC(model, 500d, 0, 0.1, 0, 0.1, 0, 0.1);
            //AddInitialConditions(model, 500d, 0, 0.1, 0, 0.1, 0, 0.1);
            
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
            dynamicAnalyzer.Solve();

            int totalNewmarkstepsNum = (int)Math.Truncate(TotalTime / TimeStep);
            var cox = new double[totalNewmarkstepsNum];
            for (int i1 = 0; i1 < totalNewmarkstepsNum; i1++)
            {
                var timeStepResultsLog = dynamicAnalyzer.ResultStorage.Logs[i1];
                cox[i1] = ((DOFSLog)timeStepResultsLog).DOFValues[model.GetNode(nodeIdToMonitor), coxMonitorDOF];
            }
            
            CSVExporter.ExportVectorToCSV(cox, "../../../Integration/Tc_nodes_mslv.csv");
        }

        private void AddTopRightBackNodesBC(Model model, double boundaryCondition,
                                                        double modelMinX,double modelMaxX,
                                                        double modelMinY,double modelMaxY,
                                                        double modelMinZ,double modelMaxZ)  
        {
            var topNodes = new List<INode>();
            var rightNodes = new List<INode>();
            var backNodes = new List<INode>();
            
            var tol = 1E-5;

            foreach (var node in model.NodesDictionary.Values)
            {
                if  (Math.Abs(modelMaxZ - node.Z) < tol) topNodes.Add(node);
                if (Math.Abs(modelMaxX - node.X) < tol) rightNodes.Add(node);
                if (Math.Abs(modelMaxY - node.Y) < tol) backNodes.Add(node);
            }
            //Union of all boundary nodes in a single enumerable.
            var peripheralNodes = topNodes.Union(backNodes).Union(rightNodes);

            var dirichletBCs = new List<NodalUnknownVariable>();

            //Add the prescribed value to all boundary nodes
            foreach (var node in peripheralNodes)
            {
                dirichletBCs.Add(new NodalUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, boundaryCondition));
            }

            model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(dirichletBCs, new INodalConvectionDiffusionNeumannBoundaryCondition[] {}));
        }
        
        private void AddInitialConditions(Model model, double initialCondition,
                                                  double modelMinX,double modelMaxX,
                                                  double modelMinY,double modelMaxY,
                                                  double modelMinZ,double modelMaxZ)
        {
            var innerBulkNodes = new List<INode>();
            var tol = 1E-5;
            var initialConditions = new List<INodalConvectionDiffusionInitialCondition>();
            foreach (var node in model.NodesDictionary.Values)
            {
                if  ((Math.Abs(modelMaxZ - node.Z) >= tol) && (Math.Abs(modelMaxX - node.X) >= tol) && (Math.Abs(modelMaxY - node.Y) >= tol))
                    initialConditions.Add(new NodalInitialUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, initialCondition));
            }
            model.InitialConditions.Add(new ConvectionDiffusionInitialConditionSet(initialConditions, new DomainInitialUnknownVariable[]{ }));
        }
    }
}
