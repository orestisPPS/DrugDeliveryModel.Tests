using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Discretization.Entities;
using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.Constitutive.Structural.Transient;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.DrugDeliveryModel.Tests.Materials;
using MGroup.FEM.Structural.Continuum;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.Solvers.Direct;

namespace MGroup.DrugDeliveryModel.Tests.EquationModels
{
	public class Eq9ModelProviderForStaggeredSolutionex8
    {
        //private double sc = 0.1;
        private double miNormal;// = 5;//KPa
        private double kappaNormal;// = 6.667; //Kpa
        private double miTumor;// = 22.44; //Kpa
        private double kappaTumor;// = 216.7; //Kpa
        private StructuralDof loadedDof;
        public double load_value { get; private set; }
    
        private double modelMinX;
        private double modelMaxX;
        private double modelMinY;
        private double modelMaxY;
        private double modelMinZ;
        private double modelMaxZ;
    
        private ComsolMeshReader reader;
       
        private Dictionary<int, double> lambda;
        Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints;

        public int nodeIdToMonitor { get; private set; } //TODO put it where it belongs (coupled7and9eqsSolution.cs)
        private StructuralDof dofTypeToMonitor;
        
        

        public int loadedNode_Id { get; private set; }

        public Eq9ModelProviderForStaggeredSolutionex8(ComsolMeshReader comsolReader, double sc, double miNormal, double kappaNormal, double miTumor,
            double kappaTumor, double timeStep, double totalTime,
            Dictionary<int, double> lambda, Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints,
            int nodeIdToMonitor, StructuralDof dofTypeToMonitor, StructuralDof loadedDof,
            double load_value, double modelMinX, double modelMaxX, double modelMinY, double modelMaxY, double modelMinZ, double modelMaxZ)
        {
            //this.sc = sc;
            this.miNormal = miNormal;
            this.kappaNormal = kappaNormal;
            this.miTumor = miTumor;
            this.kappaTumor = kappaTumor;
            this.pressureTensorDivergenceAtElementGaussPoints = pressureTensorDivergenceAtElementGaussPoints;
            this.lambda = lambda;
            reader = comsolReader;

            this.modelMinX = modelMinX;
            this.modelMaxX = modelMaxX;
            this.modelMinY = modelMinY;
            this.modelMaxY = modelMaxY;
            this.modelMinZ = modelMinZ;
            this.modelMaxZ = modelMaxZ;

            //log
            this.nodeIdToMonitor = nodeIdToMonitor;
            this.dofTypeToMonitor = dofTypeToMonitor;

            //Load
            this.loadedDof = loadedDof;
            this.load_value = load_value;

        }

        public Model GetModel()
        {
            var nodes = reader.NodesDictionary;
            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(id: 0);

            foreach (var node in nodes.Values)
            {
                model.NodesDictionary.Add(node.ID, node);
            }
            
            var materialNormal = new NeoHookeanMaterial3dJ3Isochoric(miNormal, kappaNormal);
            var materialTumor = new NeoHookeanMaterial3dJ3Isochoric(miTumor, kappaTumor);

            var elasticMaterial = new ElasticMaterial3D(youngModulus: 1, poissonRatio: 0.3);
            var DynamicMaterial = new TransientAnalysisProperties(density: 1, rayleighCoeffMass: 0, rayleighCoeffStiffness: 0);
            var elementFactory = new ContinuumElement3DFactory(elasticMaterial, DynamicMaterial);

            //var domains = new Dictionary<int, double[]>(2);
            foreach (var elementConnectivity in reader.ElementConnectivity)
            {
                var domainId = elementConnectivity.Value.Item3;
                var element = elementFactory.CreateNonLinearElementGrowt(elementConnectivity.Value.Item1, elementConnectivity.Value.Item2, domainId == 0 ? materialTumor : materialNormal, DynamicMaterial, lambda[elementConnectivity.Key]);
                element.volumeForce = pressureTensorDivergenceAtElementGaussPoints[elementConnectivity.Key][0];
                model.ElementsDictionary.Add(elementConnectivity.Key, element);
                model.SubdomainsDictionary[0].Elements.Add(element);
            }
            return model;
        }
        
        public void AddEq9ModelLoads(Model model)
        {
            INode maxDistanceNode = null;
            double currentMaxDistance = 0;
            foreach (INode node in model.NodesDictionary.Values)
            {
                double distance = Math.Sqrt(Math.Pow(node.X, 2) + Math.Pow(node.Y, 2) + Math.Pow(node.Z, 2));
                if (distance > currentMaxDistance)
                {
                    currentMaxDistance = distance;
                    maxDistanceNode = node;
                }
            }


            /*loadedNode_Id = maxDistanceNode.ID;
            nodeIdToMonitor = loadedNode_Id;*/
            var loads = new List<INodalLoadBoundaryCondition>();

            loads.Add(new NodalLoad
            (
                maxDistanceNode,
                loadedDof,
                amount: load_value
            ));

            var emptyConstraints = new List<INodalDisplacementBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(emptyConstraints, loads));
        }

        public void AddEq9ModelLoadsCenter(Model model)
        {
            var loads = new List<INodalLoadBoundaryCondition>();

            loads.Add(new NodalLoad
            (
                model.NodesDictionary[nodeIdToMonitor],
                loadedDof,
                amount: load_value
            ));

            var emptyConstraints = new List<INodalDisplacementBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(emptyConstraints, loads));

            
        }



        public void AddBottomLeftRightFrontBackBCs(Model model)
        {
            var bottomNodes = new List<INode>();
            var leftNodes = new List<INode>();
            var rightNodes = new List<INode>();
            var frontNodes = new List<INode>();
            var backNodes = new List<INode>();

            var innerBulkNodes = new List<INode>();
            var tol = 1E-5;
            //Search for all boundary nodes
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMinZ - node.Z) < tol) bottomNodes.Add(node);
                if(Math.Abs(modelMinX - node.X) < tol) leftNodes.Add(node);
                if(Math.Abs(modelMaxX - node.X) < tol) rightNodes.Add(node);
                if(Math.Abs(modelMinY - node.Y) < tol) frontNodes.Add(node);
                if(Math.Abs(modelMaxY - node.Y) < tol) backNodes.Add(node);
            }

            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMinZ - node.Z) < tol) { }

                else if (Math.Abs(modelMinX - node.X) < tol) { } 
                else if (Math.Abs(modelMaxX - node.X) < tol) { }

                else if (Math.Abs(modelMinY - node.Y) < tol) { }
                else if (Math.Abs(modelMaxY - node.Y) < tol) { }
                else innerBulkNodes.Add(node);
            }

            var constraints = new List<INodalDisplacementBoundaryCondition>();
            
            //Apply roller constraint to bottom nodes. (constrained movement in z direction)
            foreach (var node in bottomNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationZ, amount: 0d));

            //Apply roller constraint to left nodes. (constrained movement in x direction)
            foreach (var node in leftNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationX, amount: 0d));
            
            //Apply roller constraint to right nodes. (constrained movement in x direction)
            foreach (var node in rightNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationX, amount: 0d));
            
            //Apply roller constraint to front nodes. (constrained movement in y direction)
            foreach (var node in frontNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationY, amount: 0d));
            
            //Apply roller constraint to back nodes. (constrained movement in y direction)
            foreach (var node in backNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationY, amount: 0d));
            
            var emptyloads = new List<INodalLoadBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, emptyloads));
        }
        
        public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog(Model model, double pseudoTimeStep, double pseudoTotalTime, int currentStep, int nIncrements)
        {
            var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
            var algebraicModel =  solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var provider = new ProblemStructural(model, algebraicModel);
            var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, provider, numIncrements: nIncrements)
            {
                ResidualTolerance = 1E-8,
                MaxIterationsPerIncrement = 100,
                NumIterationsForMatrixRebuild = 1
            };
            var nlAnalyzer = loadControlAnalyzerBuilder.Build();
            var loadControlAnalyzer = (LoadControlAnalyzer)nlAnalyzer;
            loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
                    

                }, algebraicModel
            );

            //var analyzer = (new PseudoTransientAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, currentStep: currentStep)).Build();
            var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime,false, currentStep: currentStep);
            analyzerBuilder.SetNewmarkParametersForConstantAcceleration();
            var analyzer = analyzerBuilder.Build();


            //Sparse tet Mesh
            var watchDofs = new[]
            {
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
                    (model.NodesDictionary[nodeIdToMonitor], StructuralDof.TranslationY),
                    (model.NodesDictionary[nodeIdToMonitor], StructuralDof.TranslationZ),
                }
            };

            loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);

            return (analyzer, solver, loadControlAnalyzer);
        }
        
    }
}
