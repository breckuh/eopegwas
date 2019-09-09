#! /usr/bin/env node

const { jtree } = require("jtree")
const { Disk } = require("jtree/products/Disk.node.js")

const grammarCode = Disk.read(__dirname + `/eopClinicalSsv.grammar`)
const grammarProgram = new jtree.GrammarProgram(grammarCode)

const simulatedSsv = grammarProgram
  ._getRootNodeTypeDefinitionNode()
  .generateSimulatedData(109)
  .join("\n")

const data = jtree.TreeNode.fromSsv(
  `patientId preeclampsia momAge bmi gestAge ethnicity1 ethnicityPercent1 ethnicity2 ethnicityPercent2 ethnicity3 ethnicityPercent3 ethnicityCategory primaGravida 3HrGlucose1HrGlucola gdm iugr abruption daughter babyWgtGram placentaWgtGram\n` +
    simulatedSsv
)

Disk.write(__dirname + "/mockData/clinical.csv", data.toCsv())
