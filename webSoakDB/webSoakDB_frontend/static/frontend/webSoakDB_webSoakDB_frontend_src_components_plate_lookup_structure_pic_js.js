/*
 * ATTENTION: The "eval" devtool has been used (maybe by default in mode: "development").
 * This devtool is neither made for production nor for readable output files.
 * It uses "eval()" calls to create a separate source file in the browser devtools.
 * If you are trying to read the output file, select a different devtool (https://webpack.js.org/configuration/devtool/)
 * or disable the default devtool with "devtool: false".
 * If you are looking for production-ready output files, see mode: "production" (https://webpack.js.org/configuration/mode/).
 */
(self["webpackChunkwebSoakDB"] = self["webpackChunkwebSoakDB"] || []).push([["webSoakDB_webSoakDB_frontend_src_components_plate_lookup_structure_pic_js"],{

/***/ "./webSoakDB/webSoakDB_frontend/src/components/plate_lookup/structure_pic.js":
/*!***********************************************************************************!*\
  !*** ./webSoakDB/webSoakDB_frontend/src/components/plate_lookup/structure_pic.js ***!
  \***********************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
eval("__webpack_require__.r(__webpack_exports__);\n/* harmony export */ __webpack_require__.d(__webpack_exports__, {\n/* harmony export */   \"default\": () => (__WEBPACK_DEFAULT_EXPORT__)\n/* harmony export */ });\n/* harmony import */ var react__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! react */ \"./node_modules/react/index.js\");\n\n\nfunction StructurePic(props) {\n  var url_smiles = props.smiles.replace('#', '%23');\n  var img_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/' + url_smiles + '/PNG'; //const img_url = 'https://cactus.nci.nih.gov/chemical/structure/' + url_smiles + '/image'\n\n  return /*#__PURE__*/react__WEBPACK_IMPORTED_MODULE_0__.createElement(\"img\", {\n    className: \"img2d\",\n    src: img_url\n  });\n}\n\n/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (StructurePic);\n\n//# sourceURL=webpack://webSoakDB/./webSoakDB/webSoakDB_frontend/src/components/plate_lookup/structure_pic.js?");

/***/ })

}]);