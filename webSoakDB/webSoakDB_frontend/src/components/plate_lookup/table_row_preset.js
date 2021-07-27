import React from 'react';
import TableRowCherryPick from './table_row_cherry_pick.js';

class TableRowPreset extends TableRowCherryPick {

    getLibraryCell(){
        const classes = this.getClasses();
        return (
            <td className={classes.show_library}>
                {(this.props.compound.library !== undefined) ?
                 this.props.compound.library : 
                 '...'}
            </td>);
    }
}

export default TableRowPreset;
