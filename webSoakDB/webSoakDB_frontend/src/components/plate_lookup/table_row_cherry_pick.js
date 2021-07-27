import React from 'react';
import StructurePic from './structure_pic.js';
import TableRowPlate from './table_row_plate.js';

class TableRowCherryPick extends TableRowPlate {
  
  getGeneralCells(){
    const classes = this.getClasses();
    return (
        <React.Fragment>
          <td className={classes.show_structure}>
            <StructurePic id={this.props.compound.id} />
          </td>
          <td className={classes.show_code} > {this.props.compound.code} </td>
          <td className={classes.show_smiles} >{this.props.compound.smiles}</td>
          
        </React.Fragment>
     )
  }

  getStatValue(option){
    return this.props.compound[option];
  }
}

export default TableRowCherryPick;
