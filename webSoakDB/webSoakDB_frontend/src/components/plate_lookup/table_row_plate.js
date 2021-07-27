import React from 'react';
import {display_options, descriptor_names} from './display_options.js';
import StructurePic from './structure_pic.js';

const stat_options = display_options.slice(5)

class TableRowPlate extends React.Component {
  
  getLibraryCell(){
      return null;
  }

  getClasses(){
    let classes = {};
    
    display_options.forEach(option => {
      classes[option[0]] = this.props.display[option[0]] ? "" : "hidden" ;
    });
    return classes;
  }

  getGeneralCells(){
    const classes = this.getClasses();
    return (
        <React.Fragment>
          <td className={classes.show_well} >{this.props.compound.well}</td>
          <td className={classes.show_structure}>
            <StructurePic id={this.props.compound.compound.id} />
          </td>
          <td className={classes.show_code} > {this.props.compound.compound.code}</td>
          <td className={classes.show_smiles} >{this.props.compound.compound.smiles}</td>
          
          <td className={classes.show_concentration}>{this.props.compound.concentration}</td>
        </React.Fragment>
     )
  }

  getStatValue(option){
    return this.props.compound.compound[option];
  }

  getStatCells(){
    const classes = this.getClasses();
    return descriptor_names.map((option, index) => {
        let value = this.getStatValue(option);
        if (value % 1 !== 0){
          value = value.toFixed(2);
        }
        
        return <td key={index} className={classes[stat_options[index][0]]}>{value}</td>
      });
  }
  render() {
    const general_cells = this.getGeneralCells();
    const stat_cells = this.getStatCells();
    const library_cell = this.getLibraryCell();
    
    return (
      <tr>
        <td>{this.props.counter}</td>
        {library_cell}
        {general_cells}
        {stat_cells}
      </tr>
    );
  }
}

export default TableRowPlate;
