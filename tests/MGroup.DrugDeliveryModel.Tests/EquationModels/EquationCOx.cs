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
    public class EquationsTests13DistributedModelBuilderCox
    {
        //---------------------------------------Equation Cox---------------------------------------

        //---------------------------------------Variables------------------------------------------
        /// <summary>
        /// The average value of the three components of the fluid velocity vector  [m/s]
        /// </summary>
        private Dictionary<int, double[]> FluidSpeed = new Dictionary<int, double[]>(); // 2.32E-4 [m/s]
        double FluidInit = 2.32;
        /// <summary>
        /// Diffusivity of oxygen [m2/s]
        /// </summary>
        private const double Dox = 1.79E-4; // [m2/s]

        /// <summary>
        /// Oxygen uptake [mol/(m3*s)]
        /// </summary>
        private const double Aox = 2.55E-2; // [mol/(m3*s)]

        /// <summary>
        /// Oxygen uptake [mol/m3]
        /// </summary>
        private const double Kox = 4.64E-3; // [mol / m3]

        /// <summary>
        /// Oxygen permeability across tumor vessel walls [m/s]
        /// </summary>
        private const double PerOx = 3.55E-4; // [m/s]

        /// <summary>
        /// Vascular Density [1/m]
        /// </summary>
        private const double Sv = 7E3; // [1/m]

        /// <summary>
        /// Initial Oxygen Concentration [mol/m3]
        /// </summary>
        private const double CInitOx = 0.0; // [mol/m3]

        /// <summary>
        /// Cancer cell density [1]
        /// </summary>
        private Dictionary<int, double> T = new Dictionary<int, double>();// 500 [cells]

        //---------------------------------------Logging----------------------------------
        /// <summary>
        /// The degree of freedom that will be monitored for equation cox
        /// </summary>
        private ConvectionDiffusionDof coxMonitorDOF = ConvectionDiffusionDof.UnknownVariable;

        /// <summary>
        /// The coordinates of the monitored node
        /// </summary>
        private double[] monitorNodeCoords = { 0.09, 0.0, 0.05 };


        //---------------------------------------Time Discretization Specs------------------------------
        private const double TotalTime = 1E-3;

        /// <summary>
        /// For increased accuracy use time-step of order 1E-5
        /// </summary>
        private const double TimeStep = 1E-5;// sec

        /// <summary>
        /// Simplified version of the production term without non-linear term
        /// </summary>

        private readonly Func<double> independentLinearSource = () => PerOx * Sv * 0.2d;

        //static double[] expectedLinSolution = new double[] {-0.0002105373125587959, -0.0003133024869569006, -0.0004141663915288726, -0.0005129763178853869, -0.0006096861342003591, -0.0007042852696518675, -0.000796775015903482, -0.0008871606186806564, -0.0009754486386300758, -0.0010616460705799, -0.0011457600495923, -0.0012277977528317, -0.0013077663667799, -0.0013856730762568, -0.0014615250607187, -0.001535329492985, -0.0016070935387757, -0.0016768243565195, -0.0017445290972526, -0.0018102149045475, -0.0018738889144529, -0.0019355582554366, -0.0019952300483297, -0.0020529114062723, -0.0021086094346591, -0.0021623312310847, -0.0022140838852908, -0.0022638744791122, -0.0023117100864235, -0.002357597773087, -0.0024015445968996, -0.0024435576075406, -0.0024836438465199, -0.0025218103471266, -0.0025580641343773, -0.0025924122249651, -0.0026248616272091, -0.0026554193410037, -0.0026840923577687, -0.0027108876603993, -0.0027358122232166, -0.0027588730119183, -0.0027800769835298, -0.0027994310863556, -0.0028169422599309, -0.0028326174349734, -0.0028464635333357, -0.0028584874679578, -0.0028686961428199, -0.0028770964528957, -0.0028836952841052, -0.0028884995132695, -0.0028915160080635, -0.0028927516269715, -0.0028922132192408, -0.0028899076248369, -0.0028858416743992, -0.0028800221891956, -0.002872455981079, -0.0028631498524429, -0.002852110596178, -0.002839344995629, -0.0028248598245511, -0.0028086618470676, -0.0027907578176274, -0.0027711544809627, -0.0027498585720472, -0.0027268768160546, -0.0027022159283171, -0.0026758826142845, -0.0026478835694834, -0.0026182254794772, -0.0025869150198252, -0.0025539588560435, -0.002519363643565, -0.0024831360277003, -0.0024452826435988, -0.0024058101162096, -0.0023647250602439, -0.0023220340801361, -0.0022777437700062, -0.0022318607136227, -0.0021843914843647, -0.0021353426451856, -0.0020847207485757, -0.0020325323365265, -0.0019787839404943, -0.0019234820813643, -0.0018666332694151, -0.0018082440042839, -0.0017483207749308, -0.0016868700596048, -0.0016238983258088, -0.0015594120302662, -0.0014934176188866, -0.0014259215267321, -0.0013569301779847, -0.0012864499859127, -0.0012144873528385, -0.0011410486701057};
        static double[] expectedNonLinSolution = new double[] {-0.0002577634612139296, -0.0003824929072089781, -0.0005041287202650932, -0.0006223680518233604, -0.0007371803097879908, -0.0008485826525510643, -0.0009566085241174244, -0.0010612972094087, -0.0011626904327479, -0.0012608313031664, -0.0013557640300061, -0.0014475338808293, -0.0015361872058396, -0.0016217714711329, -0.0017043352825885, -0.0017839283953091, -0.0018606017078018, -0.0019344072414479, -0.0020053981062035, -0.0020736284535778, -0.0021391534179488, -0.0022020290472786, -0.0022623122242912, -0.0023200605791849, -0.0023753323949599, -0.0024281865064569, -0.0024786821942082, -0.0025268790742091, -0.0025728369847168, -0.0026166158711677, -0.0026582756702897, -0.0026978761944507, -0.0027354770172447, -0.0027711373612651, -0.0028049159889523, -0.0028368710973323, -0.0028670602173842, -0.00289554011869, -0.0029223667199284, -0.0029475950056838, -0.0029712789499459, -0.00299347144658, -0.0030142242469578, -0.003033587904849, -0.0030516117285885, -0.0030683437404591, -0.003083830643154, -0.0030981177931242, -0.0031112491805564, -0.0031232674156824, -0.0031342137210804, -0.0031441279295981, -0.0031530484875047, -0.0031610124624635, -0.0031680555559069, -0.0031742121193948, -0.0031795151745385, -0.0031839964360835, -0.0031876863377504, -0.0031906140604573, -0.0031928075625557, -0.003194293611741, -0.0031950978183135, -0.0031952446694936, -0.003194757564515, -0.0031936588502453, -0.0031919698571059, -0.0031897109350848, -0.0031869014896599, -0.0031835600174704, -0.0031797041415931, -0.003175350646301, -0.0031705155111971, -0.0031652139446341, -0.0031594604163457, -0.0031532686892265, -0.0031466518502135, -0.0031396223402308, -0.0031321919831694, -0.0031243720138842, -0.0031161731051972, -0.003107605393903, -0.003098678505779, -0.0030894015796072, -0.0030797832902208, -0.0030698318705901, -0.0030595551329673, -0.0030489604891118, -0.0030380549696194, -0.0030268452423817, -0.0030153376302006, -0.0030035381275875, -0.002991452416775, -0.0029790858829683, -0.0029664436288675, -0.0029535304884879, -0.0029403510403071, -0.0029269096197679, -0.002913210331164, -0.0028992570589357};
            
        public void EquationsTests13DistributedModelBuilder()
        {
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }

        private Dictionary<int, Func<double, double>> ProductionFuncsWithoutConstantTerm = new Dictionary<int, Func<double, double>>();
        public Func<double, double> getProductionFuncWithoutConstantTerm(int i)
        {
            return (double Cox) => -PerOx * Sv * Cox - Aox * T[i] * Cox / (Cox + Kox);
        }

        private Dictionary<int, Func<double, double>> ProductionFuncsWithoutConstantTermDerivative = new Dictionary<int, Func<double, double>>();
        public Func<double, double> getProductionFuncWithoutConstantTermDerivative(int i)
        {
            return (double Cox) => 0 -PerOx * Sv - Aox * T[i] / (Cox + Kox) + Aox * T[i] * Cox * Math.Pow(Cox + Kox, -2);
        }

        [Theory]
        [InlineData("../../../DataFiles/workingTetMesh2185_1Domain.mphtxt")]
        public void SolveEquationCOxNonLinearProduction(string fileName)
        {
            var mesh = new ComsolMeshReader(fileName);

            foreach (var elem in mesh.ElementConnectivity)
            {
                FluidSpeed.Add(elem.Key, new double[] { FluidInit, FluidInit, FluidInit });
                T.Add(elem.Key, 500);
                ProductionFuncsWithoutConstantTerm.Add(elem.Key, getProductionFuncWithoutConstantTerm(elem.Key));
                ProductionFuncsWithoutConstantTermDerivative.Add(elem.Key, getProductionFuncWithoutConstantTermDerivative(elem.Key));
            }
            
            var convectionDiffusionDirichletBC = new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>()
            {
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.TopDirichlet, new ConvectionDiffusionDof[] { ConvectionDiffusionDof.UnknownVariable }, new double[1][]{new double[3] {0.1, 0.1, 0.1}}, new double[] {0.2}),
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.RightDirichlet,new ConvectionDiffusionDof[] { ConvectionDiffusionDof.UnknownVariable }, new double[1][]{new double[3] {0.1, 0.1, 0.1}}, new double[] {0.2}),
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.BackDirichlet, new ConvectionDiffusionDof[] { ConvectionDiffusionDof.UnknownVariable }, new double[1][]{new double[3] {0.1, 0.1, 0.1}}, new double[] {0.2}),
            };
            var convectionDiffusionNeumannBC = new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>();

            var nodeIdToMonitor = Utilities.FindNodeIdFromNodalCoordinates(mesh.NodesDictionary, monitorNodeCoords, 1e-2);

            var modelBuilder = new CoxModelBuilder(mesh, FluidSpeed, Dox, Aox, Kox, PerOx, Sv, CInitOx, T, CInitOx, independentLinearSource, ProductionFuncsWithoutConstantTerm, ProductionFuncsWithoutConstantTermDerivative, nodeIdToMonitor, coxMonitorDOF, convectionDiffusionDirichletBC, convectionDiffusionNeumannBC);
            var model = modelBuilder.GetModel();
            modelBuilder.AddBoundaryConditions(model);

            var analysisData =  modelBuilder.GetAppropriateSolverAnalyzerAndLog(model, TimeStep, TotalTime, 0, 0);
            var dynamicAnalyzer = (NewmarkDynamicAnalyzer) analysisData.analyzer;
            dynamicAnalyzer.Initialize();
            Console.WriteLine("Solving Cox Non-Linear prod");
            dynamicAnalyzer.Solve();

            int totalNewmarkstepsNum = (int)Math.Truncate(TotalTime / TimeStep);
            var cox = new double[totalNewmarkstepsNum];
            for (int i1 = 0; i1 < totalNewmarkstepsNum; i1++)
            {
                var timeStepResultsLog = dynamicAnalyzer.ResultStorage.Logs[i1];
                cox[i1] = ((DOFSLog)timeStepResultsLog).DOFValues[model.GetNode(nodeIdToMonitor), coxMonitorDOF];
            }

            CSVExporter.ExportVectorToCSV(cox, "../../../Integration/cox_non_linear_nodes_mslv.csv");
            Console.WriteLine("FINISHED solving Cox Non-Linear prod");
            Assert.True(CompareResults(cox));

        }

        private bool CompareResults(double[] solution)
        {
            bool ret = true;
            double tolerance = 1e-1;

            for (var i=0 ; i < Math.Min(expectedNonLinSolution.Length, solution.Length) ; i++)
            {
                var error = Math.Abs(solution[i] - expectedNonLinSolution[i]);
                if ( error > tolerance)
                {
                    Console.WriteLine("Wrong result on step "+i+", Error: "+error);
                    ret = false;
                    break;
                }
            }

            return ret;
        }


    }
}
