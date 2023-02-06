using System;
using System.Collections.Generic;
using System.Linq;
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
using MGroup.FEM.ConvectionDiffusion.Tests.Commons;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Solution;
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using MGroup.Constitutive.ConvectionDiffusion.InitialConditions;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class Eq78ModelProviderForStaggeredSolutionex32
    {

        private double modelMinX;
        private double modelMaxX;
        private double modelMinY;
        private double modelMaxY;
        private double modelMinZ;
        private double modelMaxZ;

        private double Sv;
        private double k_th;
        private double pv;
        private double pl;
        private double Lp;
        private double LplSvl;

        private double boundaryValue;
        private double initialCondition;
        
        public Dictionary<int, double[]> div_vs { get; set; }
        private int nodeIdToMonitor; //TODO put it where it belongs (coupled7and9eqsSolution.cs)
        private ConvectionDiffusionDof dofTypeToMonitor;
        private ComsolMeshReader modelReader;



        public Eq78ModelProviderForStaggeredSolutionex32(ComsolMeshReader modelReader,
            double kth, double Lp, double Sv, double pv, double LplSvl, double pl, Dictionary<int, double[]> div_vs,
            double boundaryValueAllBoundaries, double initialCondition,
            int nodeIdToMonitor, ConvectionDiffusionDof dofTypeToMonitor,
            double modelMinX, double modelMaxX, double modelMinY, double modelMaxY, double modelMinZ, double modelMaxZ)
        {
            this.Sv = Sv;
            this.k_th = kth;
            this.pv = pv;
            this.pl = pl;
            this.Lp = Lp;
            this.LplSvl = LplSvl;
            this.div_vs = div_vs;

            this.modelReader = modelReader;
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;

            this.boundaryValue = boundaryValueAllBoundaries;
            this.initialCondition = initialCondition;

            this.modelMinX = modelMinX;
            this.modelMaxX = modelMaxX;
            this.modelMinY = modelMinY;
            this.modelMaxY = modelMaxY;
            this.modelMinZ = modelMinZ;
            this.modelMaxZ = modelMaxZ;

            //log
            this.nodeIdToMonitor = nodeIdToMonitor;
            this.dofTypeToMonitor = dofTypeToMonitor;
        }

        public Model GetModel()
        {
            double capacity = 0;
            double convectionCoeff = 0;
            double diffusion = k_th;
            double dependentProductionCoeff = -(Lp * Sv + LplSvl);

            Dictionary<int, double[]>ConvectionCoeffs = new Dictionary<int, double[]>(); //=> new[]  {1d, 1d, 1d};               
            Dictionary<int, double> DependentProductionCoeffs = new Dictionary<int, double>();
            Dictionary<int, double> IndependentProductionCoeffs = new Dictionary<int, double>();
            foreach (var elementConnectivity in modelReader.ElementConnectivity)
            {
                ConvectionCoeffs[elementConnectivity.Key] = new double[]
                    { convectionCoeff, convectionCoeff, convectionCoeff };
                DependentProductionCoeffs[elementConnectivity.Key] = dependentProductionCoeff;
                
                var nodes = elementConnectivity.Value.Item2;

                //var independentSource = Lp * Sv * pv  - div_vs[elementConnectivity.Key][0]; 
                //var independentSource = Lp * Sv * pv; 
                var independentSource = Lp * Sv * pv + LplSvl * pl - div_vs[elementConnectivity.Key][0]; //TODO [0] is the first gauss point Make it more genreal for all guss paints
                IndependentProductionCoeffs[elementConnectivity.Key] = independentSource;
            }

            //initialize mpdel provider solution
            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(modelReader);

            var model = modelProvider.CreateModelFromComsolFile(ConvectionCoeffs, diffusion,
                DependentProductionCoeffs, IndependentProductionCoeffs, capacity);
            return model;
        }

        public void AddAllBoundaryNodesBC(Model model)  
        {
            var topNodes = new List<INode>();
            var bottomNodes = new List<INode>();
            var leftNodes = new List<INode>();
            var rightNodes = new List<INode>();
            var frontNodes = new List<INode>();
            var backNodes = new List<INode>();
            var internalNodes = new List<INode>();
            var tol = 1E-5;
            //Search for all boundary nodes
            foreach (var node in model.NodesDictionary.Values)
            {
                if  (Math.Abs(modelMaxZ - node.Z) < tol) topNodes.Add(node);
                else if (Math.Abs(modelMinZ - node.Z) < tol) bottomNodes.Add(node);

                else if (Math.Abs(modelMinX - node.X) < tol) leftNodes.Add(node);
                else if (Math.Abs(modelMaxX - node.X) < tol) rightNodes.Add(node);

                else if (Math.Abs(modelMinY - node.Y) < tol) frontNodes.Add(node);
                else if (Math.Abs(modelMaxY - node.Y) < tol) backNodes.Add(node);
                else internalNodes.Add(node);
            }
            //Union of all boundary nodes in a single enumerable.
            var totalPeripheralNodes = leftNodes.Union(rightNodes).Union(frontNodes).Union(backNodes).Union(topNodes).Union(bottomNodes);

            var dirichletBCs = new List<NodalUnknownVariable>();

            //Add the prescribed value to all boundary nodes
            foreach (var node in totalPeripheralNodes)
            {
                dirichletBCs.Add(new NodalUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, boundaryValue));
            }

            model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(dirichletBCs, new INodalConvectionDiffusionNeumannBoundaryCondition[] {}));
        }
        
        public void AddEq78ModelInitialConditions(Model model)
        {
            var topNodes = new List<INode>();
            var bottomNodes = new List<INode>();
            var leftNodes = new List<INode>();
            var rightNodes = new List<INode>();
            var frontNodes = new List<INode>();
            var backNodes = new List<INode>();
            var innerBulkNodes = new List<INode>();
            var tol = 1E-5;
            
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMaxZ - node.Z) < tol) topNodes.Add(node);
                else if (Math.Abs(modelMinZ - node.Z) < tol) bottomNodes.Add(node);

                else if (Math.Abs(modelMinX - node.X) < tol) leftNodes.Add(node);
                else if (Math.Abs(modelMaxX - node.X) < tol) rightNodes.Add(node);

                else if (Math.Abs(modelMinY - node.Y) < tol) frontNodes.Add(node);
                else if (Math.Abs(modelMaxY - node.Y) < tol) backNodes.Add(node);
                else innerBulkNodes.Add(node);
            }

            var initialConditions = new List<INodalConvectionDiffusionInitialCondition>();
            foreach (var node in innerBulkNodes)
            {
                initialConditions.Add(new NodalInitialUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, initialCondition));
            }
            model.InitialConditions.Add(new ConvectionDiffusionInitialConditionSet(initialConditions, new DomainInitialUnknownVariable[]{ }));
        }

        

        //public void AddEq78ModelAppropriateBCs(Model model)
        //{
        //    var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(modelReader);
        //    modelProvider.AddTopAndBottomBCs(model, modelMaxZ, topValueprescribed, modelMinZ, bottomValueprescribed);
        //}

        public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog
        (Model model, double pseudoTimeStep, double pseudoTotalTime, int currentStep, int nIncrements)
        {
            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false }; //Dense Matrix Solver solves with zero matrices!
            //var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var provider = new ProblemConvectionDiffusion(model, algebraicModel);

            var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, provider, numIncrements: 1)
            {
                ResidualTolerance = 1E-8,
                MaxIterationsPerIncrement = 100,
                NumIterationsForMatrixRebuild = 1
            };
            var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();

            loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(new List<(INode node, IDofType dof)>()
            {(model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor)}, algebraicModel);

            var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, false, currentStep: currentStep);
            analyzerBuilder.SetNewmarkParametersForConstantAcceleration();
            var analyzer = analyzerBuilder.Build();
            var watchDofs = new[]
            {
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
                }
            };
            loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);

            return (analyzer, solver, loadControlAnalyzer);
        }
    }
}
