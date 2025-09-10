#!/usr/bin/env python3

import json
import argparse
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, Image, ListFlowable, ListItem
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

# Initialize styles
styles = getSampleStyleSheet()
styles.add(ParagraphStyle(
    name='CustomTitle',
    parent=styles['Title'],
    fontSize=24,
    spaceAfter=30
))
styles.add(ParagraphStyle(
    name='WarningHeader',
    parent=styles['Heading2'],
    fontSize=14,
    spaceAfter=12,
    textColor=colors.red
))
styles.add(ParagraphStyle(
    name='BoldText',
    parent=styles['BodyText'],
    fontName='Helvetica-Bold'
))
styles.add(ParagraphStyle(
    name='Indent',
    parent=styles['BodyText'],
    leftIndent=24
))

def create_table_data(data, title):
    """Convert dictionary data to table format with title"""
    table_data = [[Paragraph(f"<b>{title}</b>", styles['Heading2']), '']]  # Title row
    table_data.extend([[key, value] for key, value in data.items()])
    return table_data

def create_styled_table(table_data):
    """Create a styled table from data"""
    table = Table(table_data, colWidths=[4*inch, 2.5*inch])
    
    style = TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('SPAN', (0, 0), (-1, 0)),
        ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
        ('BACKGROUND', (0, 1), (-1, -1), colors.white),
        ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),
        ('ALIGN', (0, 1), (0, -1), 'LEFT'),
        ('ALIGN', (1, 1), (1, -1), 'RIGHT'),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 14),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 1), (-1, -1), 10),
        ('TOPPADDING', (0, 1), (-1, -1), 6),
        ('BOTTOMPADDING', (0, 1), (-1, -1), 6),
    ])
    table.setStyle(style)
    return table

def create_text_section(title, content):
    """Create flowables for a text section"""
    section = []
    
    # Add section header
    if title:
        section.append(Paragraph(title, styles['WarningHeader']))
    
    # Process content items
    for item_key, item_value in content.items():
        # Add main item in bold
        section.append(Paragraph(f"<b>{item_key}</b>", styles['BoldText']))
        
        # Handle different value types
        if isinstance(item_value, list):
            # Create bullet list for multiple values
            bullet_items = []
            for sub_item in item_value:
                bullet_items.append(ListItem(Paragraph(sub_item, styles['BodyText'])))
            section.append(ListFlowable(bullet_items, bulletType='bullet'))
        elif isinstance(item_value, dict):
            # Handle nested dictionaries
            for sub_key, sub_value in item_value.items():
                section.append(Paragraph(f"• {sub_key}: {sub_value}", styles['Indent']))
        else:
            # Handle single string values
            section.append(Paragraph(item_value, styles['Indent']))
        
        section.append(Spacer(1, 8))
    
    section.append(Spacer(1, 24))
    return section

def create_report(json_data, output_file, image_paths=None):
    doc = SimpleDocTemplate(
        output_file,
        pagesize=letter,
        rightMargin=72,
        leftMargin=72,
        topMargin=72,
        bottomMargin=72
    )
    
    story = []
    title = Paragraph("Genome Analysis Report", styles['Title'])
    story.append(title)
    story.append(Spacer(1, 30))
    
    # Process each section in JSON data
    for section_title, section_content in json_data.items():
        # Handle text sections
        if section_title.lower() == "warnings":
            story.extend(create_text_section(section_title, section_content))
        # Handle other sections as tables
        else:
            table_data = create_table_data(section_content, section_title)
            table = create_styled_table(table_data)
            story.append(table)
            story.append(Spacer(1, 30))
    
    # Add images
    if image_paths:
        for img_path in image_paths:
            try:
                img = Image(img_path, width=6*inch, height=4*inch)
                story.append(img)
                story.append(Spacer(1, 20))
            except Exception as e:
                print(f"Error adding image {img_path}: {e}")
    
    doc.build(story)

def main():
    parser = argparse.ArgumentParser(description='Convert JSON to PDF report')
    parser.add_argument('-j', '--json', required=True, help='Path to input JSON file')
    parser.add_argument('-i', '--images', help='Text file with image paths')
    parser.add_argument('-o', '--output', default='AnnoAudit_Report.pdf', help='Output PDF file path')
    
    args = parser.parse_args()

    with open(args.json) as f:
        json_data = json.load(f)
    
    image_paths = []
    if args.images:
        with open(args.images, 'r') as f:
            image_paths = [line.strip() for line in f if line.strip()]
    
    create_report(json_data, args.output, image_paths)
    print(f"Report generated: {args.output}")

if __name__ == "__main__":
    main()